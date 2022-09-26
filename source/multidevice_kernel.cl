
/*******************************************************************************************************************************************
* Matrix free projectors for various ML and MAP algorithms. This function calculates the sensitivity image d_Summ = sum(A,1) (sum of every 
* row) and rhs = A*(y./(A'*x)) (forward and backward projections), where A is the system matrix, y the measurements and x the 
* estimate/image.
* 
* This code is also used for the precomputation phase, where the number of voxels that each LOR traverses is computed.
*
* Used by implementations 2 and 3.
*
* This file contains all the three different projectors (Siddon, orthogonal, volume-based). Preprocessor commands are used to separate
* different areas of code for the different projectors. This code also includes the precomputation phase where the number of voxels in each 
* LOR are computed. Furthermore the forward-backward projection example uses this same file. 64-bit atomics are also currently included in 
* the same file and used if supported.
*
* Compiler preprocessing is utilized heavily, for example all the corrections are implemented as compiler preprocesses. The code for 
* specific correction is thus only applied if it has been selected. The kernels are always compiled on-the-fly, though when using same input 
* parameters the kernel should be loaded from cache leading to a slightly faster startup time.
*
* INPUTS:
* global_factor = a global correction factor, e.g. dead time
* d_epps = a small constant to prevent division by zero,
* d_N = d_Nx * d_Ny * d_Nz,
* d_Nx/y/z = image size in x/y/z- dimension,
* d_dz/x/y = distance between adjecent voxels in z/x/y-dimension,
* d_bz/x/y = distance from the pixel space to origin (z/x/y-dimension),
* d_bzb = part in parenthesis of equation (9) in [1] precalculated when k = Nz,
* d_maxxx/yy = maximum distance of the pixel space from origin in x/y-dimension,
* d_zmax = maximum value of d_zdet,
* d_NSlices = the number of image slices,
* d_size_x = the number of detector elements,
* d_TotSinos = Total number of sinograms,
* d_det_per_ring = number of detectors per ring,
* d_raw = if 1 then raw data is used otherwise sinogram data
* pRows = number of pseudo rings,
* d_Nxy = d_Nx * d_Ny,
* fp = if 1, then only forward projection is computed, if 2 only backprojection, if 0 then both,
* tube_width = the width of of the strip used for orthogonal distance based projector (2D),
* crystal_size_z = the width of of the tube used for orthogonal distance based projector (3D),
* bmin = smaller orthogonal distances than this are fully inside the TOR, volume projector only,
* bmax = Distances greater than this do not touch the TOR, volume projector only,
* Vmax = Full volume of the spherical "voxel", volume projector only,
* d_epsilon_mramla = Epsilon value for MRAMLA and MBSREM,
* d_TOFCenter = Offset of the TOF center from the first center of the FOV,
* d_atten = attenuation data (images),
* d_norm = normalization coefficients,
* d_scat = scatter coefficients when using the system matrix method, 
* d_Summ = buffer for d_Summ,
* d_lor = number of pixels that each LOR traverses,
* d_pseudos = location of pseudo rings,
* d_x/y/z_det = detector x/y/z-coordinates,
* x/y/z_center = Cartesian coordinates for the center of the voxels (x/y/z-axis),
* d_xy/zindex = for sinogram format they determine the detector indices corresponding to each sinogram bin (unused with raw data),
* V = Precomputed volumes for specific orthogonal distances, volume projector only
* d_L = detector numbers for raw data (unused for sinogram format),
* d_Sino = Sinogram/raw data,
* d_sc_ra = Randoms and/or scatter data,
* d_OSEM = buffer for OSEM/MLEM estimates,
* d_output = buffer for OSEM/MLEM RHS elements,
* no_norm = If 1, normalization constant is not computed,
* m_size = Total number of LORs for this subset,
* cumsum = offset for input vector b in backprojection
*
* OUTPUTS:
* d_output = RHS values for OSEM/MLEM,
* d_OSEM = OSEM/MLEM estimate,
* d_Summ = Sensitivity image
*
* [1] Jacobs, F., Sundermann, E., De Sutter, B., Christiaens, M. Lemahieu, I. (1998). A Fast Algorithm to Calculate the Exact Radiological 
* Path through a Pixel or Voxel Space. Journal of computing and information technology, 6 (1), 89-94.
*
* Copyright (C) 2019-2022 Ville-Veikko Wettenhovi
*
* This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it wiL be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************************************************************************/

// Matrix free orthogonal distance-based ray tracer, no precomputation step
#if (defined(CT) || defined(PET) || defined(SPECT)) && !defined(LISTMODE)
__kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, 1)))
#else
__kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(LOCAL_SIZE, 1, 1)))
#endif
#ifdef FIND_LORS
void siddon_precomp(const uint d_Nxy, const uint d_N, const uint3 d_Nxyz, const float3 d_d, const float3 b, const float3 d_max,
	const uint d_size_x, const uint d_det_per_ring, const __global float* restrict d_xy, 
	const __global float* restrict d_z, __global ushort* restrict d_lor, const __global ushort* restrict d_L, const ulong m_size) {

#else

void kernel_multi(const float global_factor, const float d_epps, const uint d_N, const uint3 d_Nxyz, const float3 d_d, 
	const float3 b, const float3 d_max, const uint d_size_x, const uint d_det_per_ring,
	const uint d_Nxy, const uchar fp, const float sigma_x, const float2 crystalSize, const float orthWidth, const float bmin,
	const float bmax, const float Vmax,	
#ifdef MASKFP
	__read_only image2d_t maskFP,
#endif
#ifdef TOF
	__constant float* TOFCenter, 
#endif
#ifdef ORTH
	__constant float* x_center, __constant float* y_center, __constant float* z_center,
	__constant float* V, 
#endif
	__constant uchar* MethodList, const float d_epsilon_mramla,
#if !defined(CT) && defined(ATN)
	__read_only image3d_t d_atten,
#endif
#if defined(CT) || defined(SPECT) || defined(PET)
	const uint d_sizey, const long d_nProjections,
#endif
#if defined(LISTMODE)
	const __global float* restrict d_xy, const __global float* restrict d_z,
#else
	__constant float* d_xy, __constant float* d_z,
#endif
	const __global float* restrict d_norm, const __global float* restrict d_scat, __global CAST* restrict d_Summ, 
#ifdef PRECOMPUTE
	const __global ushort* restrict d_lor, 
#endif
#ifdef SUBSETS
	const __global uint* restrict d_xyindex, const __global ushort* restrict d_zindex, 
#endif
#ifdef RAW
	const __global ushort* restrict d_L, 
#endif
	const __global float* restrict d_Sino, const __global float* restrict d_sc_ra, 
#if defined(FP)
	__read_only image3d_t d_OSEM, 
#else
	const __global float* restrict d_OSEM, 
#endif
	//////////////////////////// Not yet implemented ////////////////////////////
#ifdef MATRIX
	__global int* restrict indices, __global int* restrict rowInd, __global float* restrict values, const int maxLOR, 
#endif
	//////////////////////////// Not yet implemented ////////////////////////////
#ifndef MBSREM
	__global CAST* d_output, const uchar no_norm, const ulong m_size, 
#else
	const uint d_alku, const uchar MBSREM_prepass, __global float* restrict d_ACOSEM_lhs, __global float* restrict d_Amin, __global CAST* restrict d_co,
	__global CAST* restrict d_aco, __global float* restrict d_E, const ulong m_size, const RecMethodsOpenCL MethodListOpenCL, 
#endif
	const ulong cumsum) {
#endif
	// Get the current global index
	//size_t idx = get_global_id(0);
	int3 i = { get_global_id(0), get_global_id(1), get_global_id(2) };
#if defined(CT) || defined(SPECT) || defined(PET)
	size_t idx = get_global_id(0) + get_global_id(1) * d_size_x + get_global_id(2) * d_sizey * d_size_x;
	if (i.x >= d_size_x || i.y >= d_sizey || i.z >= d_nProjections)
#else
	//size_t idx = get_global_id(0) + get_global_id(1) * get_global_size(0) + get_global_id(2) * get_global_size(1) * get_global_size(0);
	size_t idx = get_global_id(0);
	//if (idx == 1 || idx == 0) {
	//	printf("idx = %d\n", idx);
	//}
	if (idx >= m_size)
#endif
		return;


#ifdef MASKFP
//#if !defined(PET) && !defined(CT) && !defined(SPECT) && !defined(LISTMODE)
//	i.x = i.y * d_size_x;
//#endif
	const int maskVal = read_imageui(maskFP, sampler_MASK, (int2)(i.x, i.y)).x;
	if (maskVal == 0)
		return;
#endif
#ifndef FIND_LORS // Not the precomputation phase
#ifdef TOF
	float local_sino = 0.f;
#ifndef LISTMODE2
#if defined(BP) && defined(FP)
#pragma unroll NBINS
	for (long to = 0L; to < NBINS; to++)
		local_sino += d_Sino[idx + m_size * to];
#endif
#endif
#else
#ifdef LISTMODE2
	const float local_sino = 0.f;
#else
#if defined(BP) && defined(FP)
	const float local_sino = (d_Sino[idx]);
#else
	const float local_sino = 0.f;
#endif
#endif
#endif
#if !defined(MBSREM) && defined(BP) && defined(FP)
	if (no_norm == 1u && local_sino == 0.f)
		return;
#elif defined(MBSREM)
	const uchar no_norm = 0u;
#endif
#endif
	//return;

	float3 s, d;
	// Load the next detector index
#if (defined(CT) || (defined(SPECT))) && !defined(LISTMODE) && !defined(PET)
	getDetectorCoordinatesCT(d_xy, d_z, &s, &d, i, d_size_x, d_sizey, crystalSize);
#elif defined(LISTMODE)
	getDetectorCoordinatesListmode(d_xy, &s, &d, idx);
#elif defined(RAW) // raw data
	getDetectorCoordinatesRaw(d_xy, d_z, d_L, d_det_per_ring, idx, &s, &d); // Sinogram data
#elif defined(FIND_LORS) || !defined(SUBSETS) // Precomputation phase
	getDetectorCoordinatesFullSinogram(d_size_x, i, &s, &d, d_xy, d_z);
#else // Not the precomputation phase
	getDetectorCoordinates(d_xyindex, d_zindex, idx, &s, &d, d_xy, d_z);
#endif

	//if (i.x == 100 && i.y == 100 && i.z == 0) {
	//	printf("i.x = %d\n", i.x);
	//	printf("i.y = %d\n", i.y);
	//	printf("i.z = %d\n", i.z);
	//	printf("d_size_x = %d\n", d_size_x);
	//	//printf("d_sizey = %d\n", d_sizey);
	//	//printf("d_nProjections = %d\n", d_nProjections);
	//	//printf("Np_n = %d\n", Np_n);
	//	//printf("tempi = %d\n", tempi);
	//	//printf("tempj = %d\n", tempj);
	//	//printf("tempk = %d\n", tempk);
	//	//printf("get_global_size(0) = %d\n", get_global_size(0));
	//	//printf("get_global_size(1) = %d\n", get_global_size(1));
	//	//printf("get_global_size(2) = %d\n", get_global_size(2));
	//	//printf("local_sino = % f\n", local_sino);
	//	printf("s.x = % f\n", s.x);
	//	printf("s.y = % f\n", s.y);
	//	printf("s.z = % f\n", s.z);
	//	printf("d.x = % f\n", d.x);
	//	printf("d.y = % f\n", d.y);
	//	printf("d.z = % f\n", d.z);
	//	//printf("axOSEM = % f\n", axOSEM);
	//}
	// Calculate the x, y and z distances of the detector pair
	float3 diff = d - s;
	//float y_diff = (yd - ys);
	//float x_diff = (xd - xs);
	//const float z_diff = (zd - zs);

#ifdef PRECOMPUTE // Using precomputed data
	uint Np = convert_uint(d_lor[idx]);
#if !defined(DEC) || defined(TOF) // Intermediate results are not saved
	uint Np_n = Np;
#endif
#else // No precomputation
	if (all(diff == 0.f) || (diff.x == 0.f && diff.y == 0.f))
		return;
	uint Np = 0u;
#if !defined(DEC) || defined(TOF) // Intermediate results are not saved
	uint Np_n = 0u;
#endif
#endif

	bool RHS = false, SUMMA = false;
	float L = 0.f;

#ifdef ATN // Attenuation included
	float jelppi = 0.f;
#endif

	float local_norm = 0.f;
	float local_scat = 0.f;

#ifdef NORM // Normalization included
	local_norm = d_norm[idx];
#endif
#ifdef SCATTER // Scatter data included
	local_scat = d_scat[idx];
#endif

	uint d_N0 = d_Nxyz.x;
	uint d_N1 = d_Nxyz.y;
	uint d_N2 = 1u;
	uint d_N3 = d_Nxyz.x;

//#ifdef ORTH // 3D Orthogonal
//	uint d_N4 = d_Nxyz.z;
//#endif

	//int3 tempijk = { 0, 0, 0 }; 
	//int3 u = { 0, 0, 0 };
	 int tempi = 0, tempj = 0, tempk = 0, ux = 0, uy = 0, uz = 0;
	//int tempijk.x = 0, tempijk.y = 0, tempijk.z = 0, ux = 0, uy = 0, uz = 0;

#ifndef FIND_LORS // Not the precomputation phase

#ifdef AF
#ifdef MBSREM
#ifdef TOF
	float axACOSEM[NBINS];
	float ax[NBINS];
#pragma unroll NBINS
	for (uint to = 0; to < NBINS; to++) {
		axACOSEM[to] = 0.f;
		ax[to] = 0.f;
	}
#else
	float axACOSEM = 0.f;
	float axCOSEM = 0.f;
#endif
#ifdef TOF
	float minimi[NBINS];
#pragma unroll NBINS
	for (uint to = 0; to < NBINS; to++)
		minimi[to] = 1e8f;
#else
	float minimi = 1e8f;
#endif
#else
	float ax[NROLLS];
#pragma unroll
	for (uint kk = 0; kk < NROLLS; kk++)
		ax[kk] = 0.f;
#if defined(BP) && !defined(FP)
#pragma unroll NROLLS
	for (uint to = 0; to < NROLLS; to++)
		ax[to] = d_OSEM[idx + to * m_size];
#endif
#endif
#else
#ifdef TOF
	float ax[NBINS];
#pragma unroll NBINS
	for (uint to = 0; to < NBINS; to++)
		ax[to] = 0.f;
#else
	float axOSEM = 0.f;
#endif
#endif

#if !defined(FP) && defined(BP)
#ifndef LISTMODE2
	//if (fp == 2) {
#ifdef TOF
#pragma unroll NBINS
		for (uint to = 0; to < NBINS; to++)
			ax[to] = d_OSEM[idx + to * m_size + cumsum];
#else
		axOSEM = d_OSEM[idx + cumsum];
#endif
	//}
#endif
#endif

#ifdef ORTH // Orthogonal or volume-based ray tracer
	//uchar xyz = 0u;
	float kerroin = 0.f;
	__constant float* xcenter = x_center;
	__constant float* ycenter = y_center;
#if !defined(VOL) // Orthogonal
	kerroin = length(diff) * orthWidth;
#elif defined VOL // Volume-based
	kerroin = length(diff);
//#elif defined(CRYSTZ) && !defined(VOL) // 3D Orthogonal
//	kerroin = length(diff) * crystalSize.y;
#endif
#elif defined SIDDON // Siddon
	uint local_ind = 0u;
	//int indX = 0, indY = 0, indZ = 0;
	//int4 testi = { &indX, &indY, &indZ, 0 };
	int4 localInd = { 0, 0, 0, 0 };
	//uint local_ind = 0;
	float local_ele = 0.f;
#ifdef TOF
	float D = 0.f, DD = 0.f;
#endif
#endif
#else // Precomputation phase
	ushort temp_koko = 0u;
#endif
	//return;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//If the LOR is perpendicular in the y-direction (Siddon cannot be used)
	if (fabs(diff.z) < 1e-6f && (fabs(diff.y) < 1e-6f || fabs(diff.x) < 1e-6f)) {

		//return;
#ifdef FIND_LORS // Precomputation phase
		if (fabs(diff.y) < 1e-6f && d.y <= d_max.y && d.y >= b.y && s.y <= d_max.y && s.y >= b.y) {
			d_lor[idx] = convert_ushort(d_Nxyz.x);
			return;
		}
		else if (fabs(diff.x) < 1e-6f && d.x <= d_max.x && d.x >= b.x && s.x <= d_max.x && s.x >= b.x) {
			d_lor[idx] = convert_ushort(d_Nxyz.y);
			return;
		}
		else
			return;
#else // Not the precomputation phase

		tempk = convert_int(fabs(s.z - b.z) / d_d.z);
		//if (i.z == 28) {
		//	printf("tempk = %d\n", tempk);
		//}

#if defined(ORTH)
		float center2;
		float temp = 0.f;
#endif
		float d_b, dd, d_db, d_d2;
		if (fabs(diff.y) < 1e-6f && d.y <= d_max.y && d.y >= b.y && s.y <= d_max.y && s.y >= b.y) {
			d_b = b.y;
			dd = d.y;
			d_db = d_d.y;
#if defined(ORTH)
			center2 = x_center[0];
			xcenter = y_center;
#elif defined(SIDDON)
			d_d2 = d_d.x;
#endif
			float xs_apu = s.x;
			s.x = s.y;
			s.y = xs_apu;
			float xdiff_apu = diff.x;
			diff.x = diff.y;
			diff.y = xdiff_apu;
			d_N0 = d_Nxyz.y;
			d_N1 = d_Nxyz.x;
			d_N2 = d_Nxyz.y;
			d_N3 = 1u;
		}
		else if (fabs(diff.x) < 1e-6f && d.x <= d_max.x && d.x >= b.x && s.x <= d_max.x && s.x >= b.x) {
			d_b = b.x;
			dd = d.x;
#if defined(ORTH)
			center2 = y_center[0];
#elif defined(SIDDON)
			d_d2 = d_d.y;
#endif
			d_db = d_d.x;
		}
		else
			return;
#ifdef SIDDON // Siddon
		float templ_ijk = 0.f;
		//printf("templ_ijk = %f\n", templ_ijk);
		localInd.z = tempk;
		//uint z_loop = 0u;
		perpendicular_elements(d_b, d_db, d_N0, dd, d_d2, d_N1, &templ_ijk, &localInd, &tempk, d_N2, d_N3, d_norm, idx, global_factor,
			d_scat
#if !defined(CT) && defined(ATN)
			, d_atten
#endif
		);
		//local_ind = tempk;
		//perpendicular_elements(d_b, d_d, d_N0, dd, d_d2, d_N1, d_atten, &templ_ijk, &z_loop, tempk, d_N2, d_N3,
		//	d_norm, idx, global_factor, d_scat);
#ifdef TOF
		float dI = (d_d2 * d_N1) / 2.f * -sign(diff.y);
		D = dI;
		DD = D;
		float temp = templ_ijk / d_d2;
#ifdef DEC
		__private float store_elements[DEC * NBINS];
#else
		__private float store_elements[1];
#endif
#if defined(DEC) || defined(FP)
		//local_ind = z_loop;
		for (uint ii = 0; ii < d_N1; ii++) {
			const float TOFSum = TOFLoop(DD, d_d2, store_elements, TOFCenter, sigma_x, &D, ii * NBINS, d_epps);
#ifdef FP 
			denominatorTOF(ax, d_d2, d_OSEM, localInd, TOFSum, store_elements, DD, TOFCenter, sigma_x, &D, ii * NBINS, d_epps, d_N);
			//denominatorTOF(ax, d_d2, d_OSEM, local_ind, TOFSum, store_elements, DD, TOFCenter, sigma_x, &D, ii * NBINS, d_epps, d_N);
#endif
			if (d_N3 == 1)
				localInd.x++;
			else
				localInd.y++;
			//local_ind += d_N3;
		}
#endif
//#ifndef AF
#if defined(FP) && !defined(BP)
		//if (fp == 1) {
			nominatorTOF(MethodList, ax, d_Sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx, m_size, local_sino);
#pragma unroll NBINS
			for (long to = 0; to < NBINS; to++)
				d_output[idx + to * m_size] = ax[to];
			return;
		//}
#endif
//#endif
#ifdef BP
#ifndef DEC
		D = DD;
#endif
		local_ele = templ_ijk;
		//localInd.z = tempk;
		local_ind = tempk;
#if defined(FP) && defined(BP)
		if (local_sino != 0.f) {
#endif
			nominatorTOF(MethodList, ax, d_Sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx, m_size, local_sino);
			for (uint ii = 0u; ii < d_N1; ii++) {
#ifndef DEC
				const float TOFSum = TOFLoop(DD, d_d2, store_elements, TOFCenter, sigma_x, &D, ii * NBINS, d_epps);
#endif
				//backprojectTOF(local_ind, local_ele, ii * NBINS, store_elements, ax, d_Summ, 
				backprojectTOF(local_ind, localInd, local_ele, ii * NBINS, store_elements, ax, d_Summ,
#ifndef DEC
					temp, sigma_x, &D, DD, TOFCenter, d_epps, TOFSum,
#endif
#ifdef MBSREM
					MethodListOpenCL, d_alku, MBSREM_prepass, minimi, axACOSEM, d_OSEM, d_E, d_co, d_aco, local_sino, idx, m_size);
#else
					d_output, no_norm, d_N);
#endif
				local_ind += d_N3;
				if (d_N3 == 1)
					localInd.x++;
				else
					localInd.y++;
			}
#if defined(FP) && defined(BP)
		}
		else {
			for (uint ii = 0u; ii < d_N1; ii++) {
#ifndef DEC
				const float TOFSum = TOFLoop(DD, d_d2, store_elements, TOFCenter, sigma_x, &D, ii * NBINS, d_epps);
#endif
				//sensTOF(local_ind, local_ele, ii * NBINS, store_elements, d_Summ, 
				sensTOF(local_ind, localInd, local_ele, ii * NBINS, store_elements, d_Summ,
#ifndef DEC
					temp, sigma_x, &D, DD, TOFCenter, d_epps, TOFSum,
#endif
#ifdef MBSREM
					MethodListOpenCL, d_alku, MBSREM_prepass, minimi, axACOSEM, d_OSEM, d_E, idx, m_size,
#endif
					no_norm);
				local_ind += d_N3;
				if (d_N3 == 1)
					localInd.x++;
				else
					localInd.y++;
			}
		}
#endif
#ifdef MBSREM
		if (d_alku == 0u && (MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1)
#pragma unroll NBINS
			for (long to = 0L; to < NBINS; to++) {
				d_Amin[idx + to * m_size] = minimi[to];
			}
		else if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u) {
#pragma unroll NBINS
			for (long to = 0L; to < NBINS; to++) {
#ifdef RANDOMS
				axACOSEM[to] += d_sc_ra[idx];
#endif
				d_ACOSEM_lhs[idx + to * m_size] = axACOSEM[to];
			}
		}
#endif
#endif
#else
#ifdef MBSREM
		if (d_alku == 0u && ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0 ||
			(MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1 || MethodListOpenCL.RBIOSL == 1 || MethodListOpenCL.RBI == 1))) && local_sino != 0.f) {
			local_ele = templ_ijk;
			//local_ind = z_loop;
			//local_ind = z_loop;
			for (uint ii = 0u; ii < Np; ii++) {
				if (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0)
					axCOSEM += (local_ele * read_imagef(d_OSEM, samplerIm, localInd).x);
					//axCOSEM += (local_ele * d_OSEM[local_ind]);
				if (MBSREM_prepass == 1)
#ifdef ATOMIC // 64-bit atomics
					atom_add(&d_Summ[local_ind], convert_long(local_ele * TH));
					//atom_add(&d_Summ[local_ind], convert_long(local_ele * TH));
#elif defined(ATOMIC32)
					//atomic_add(&d_Summ[local_ind], convert_int(local_ele * TH));
					atomic_add(&d_Summ[local_ind], convert_int(local_ele * TH));
#else // 32-bit float atomics
					atomicAdd_g_f(&d_Summ[local_ind], local_ele);
					//atomicAdd_g_f(&d_Summ[local_ind], local_ele);
#endif
				if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
					if (local_ele < minimi && local_ele > 0.f)
						minimi = local_ele;
					d_E[idx] += local_ele;
				}
				local_ind += d_N3;
				//local_ind += d_N3;
			}
			if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1)
				d_Amin[idx] = minimi;
			if (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) {
#ifndef CT
				if (axCOSEM < d_epps)
					axCOSEM = d_epps;
#endif
#ifdef RANDOMS
				axCOSEM += d_sc_ra[idx];
#endif
#ifdef CT
				axCOSEM = native_exp(-axCOSEM) / local_sino;
#else
				axCOSEM = local_sino / axCOSEM;
#endif
			}
			local_ele = templ_ijk;
			//local_ind = z_loop;
			local_ind = tempk;
			for (uint ii = 0u; ii < Np; ii++) {
				if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2))
#ifdef ATOMIC // 64-bit atomics
					atom_add(&d_co[local_ind], convert_long(axCOSEM * local_ele * TH));
					//atom_add(&d_co[local_ind], convert_long(axCOSEM * local_ele * TH));
#elif defined(ATOMIC32)
					atomic_add(&d_co[local_ind], convert_int(axCOSEM * local_ele * TH));
					//atomic_add(&d_co[local_ind], convert_int(axCOSEM * local_ele * TH));
#else // 32-bit float atomics
					atomicAdd_g_f(&d_co[local_ind], axCOSEM * local_ele);
					//atomicAdd_g_f(&d_co[local_ind], axCOSEM * local_ele);
#endif
				if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1))
#ifdef ATOMIC // 64-bit atomics
					atom_add(&d_aco[local_ind], convert_long(axCOSEM * local_ele * TH));
					//atom_add(&d_aco[local_ind], convert_long(axCOSEM * local_ele * TH));
#elif defined(ATOMIC32)
					atomic_add(&d_aco[local_ind], convert_int(axCOSEM * local_ele * TH));
					//atomic_add(&d_aco[local_ind], convert_int(axCOSEM * local_ele * TH));
#else // 32-bit float atomics
					atomicAdd_g_f(&d_aco[local_ind], axCOSEM * local_ele);
					//atomicAdd_g_f(&d_aco[local_ind], axCOSEM * local_ele);
#endif
				local_ind += d_N3;
				//local_ind += d_N3;
			}
		}
		else if (MBSREM_prepass == 1 && d_alku == 0u) {
			local_ele = templ_ijk;
			local_ind = tempk;
			//local_ind = z_loop;
			for (uint ii = 0u; ii < Np; ii++) {
#ifdef ATOMIC // 64-bit atomics
				atom_add(&d_Summ[local_ind], convert_long(local_ele * TH));
				//atom_add(&d_Summ[local_ind], convert_long(local_ele * TH));
#elif defined(ATOMIC32)
				atomic_add(&d_Summ[local_ind], convert_int(local_ele* TH));
				//atomic_add(&d_Summ[local_ind], convert_int(local_ele* TH));
#else // 32-bit float atomics
				atomicAdd_g_f(&d_Summ[local_ind], local_ele);
				//atomicAdd_g_f(&d_Summ[local_ind], local_ele);
#endif
				if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1)) {
					if (local_ele < minimi && local_ele != 0.f)
						minimi = local_ele;
					d_E[idx] += local_ele;
				}
				local_ind += d_N3;
			}
			if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1))
				d_Amin[idx] = minimi;
		}
		else if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u) {
			local_ele = templ_ijk;
			localInd.z = tempk;
			//local_ind = z_loop;
			for (uint ii = 0u; ii < Np; ii++) {
				axACOSEM += (local_ele * read_imagef(d_OSEM, samplerIm, localInd).x);
				//axACOSEM += (local_ele * d_OSEM[local_ind]);
				if (d_N3 == 1)
					localInd.x++;
				else
					localInd.y++;
				//local_ind += d_N3;
			}
#ifdef RANDOMS
				axACOSEM += d_sc_ra[idx];
#endif
#ifdef CT
			d_ACOSEM_lhs[idx] = native_exp(-axACOSEM);
#else
			d_ACOSEM_lhs[idx] = axACOSEM;
#endif
		}
#else
#if !defined(BP) && defined(FP)
		//if (fp == 1) {
			//L = length(diff);
			//local_ele = templ_ijk / L;
			local_ele = templ_ijk;
			//localInd.z = tempk;
			//local_ind = z_loop;
			//if (i.x == 1 && i.y == 300 && i.z == 59) {
			//}
			for (uint ii = 0u; ii < d_N1; ii++) {

#ifdef AF // Implementation 2

				denominator(local_ele, ax, localInd, d_N, d_OSEM);
				//denominator(local_ele, ax, local_ind, d_N, d_OSEM);

#else // Implementation 3
				//denominator_multi(local_ele, &axOSEM, &d_OSEM[local_ind]);
				denominator_multi(local_ele, &axOSEM, read_imagef(d_OSEM, samplerIm, localInd).x);
#endif
				if (d_N3 == 1)
					localInd.x++;
				else
					localInd.y++;
				//if (d_N3 == 1)
				//	indX++;
				//else
				//	indY++;
				//	localInd.x = localInd.x + 1;
					//localInd.y = localInd.y + 1;
				//local_ind += d_N3;
					//localInd++;
			}
#ifdef AF // Implementation 2

			nominator(MethodList, ax, local_sino, d_epsilon_mramla, d_epps, 1.f, d_sc_ra, idx); 
			forwardProjectAF(d_output, ax, idx, d_N);

#else // Implementation 3
			nominator_multi(&axOSEM, local_sino, d_epps, 1.f, d_sc_ra, idx);
			d_output[idx] = axOSEM;
#endif
			return;
		//}
#endif
		//return;
		//local_ele = templ_ijk;
		//localInd.z = tempk;
		//local_ind = z_loop;
		 //Calculate the next index and store it as weL as the probability of emission
		 //If measurements are present, calculate the 
#if defined(FP) && defined(BP) // Forward projection
		if (local_sino != 0.f) {

			for (uint ii = 0u; ii < d_N1; ii++) {

#ifdef AF // Implementation 2

				denominator(local_ele, ax, localInd, d_N, d_OSEM);
				//denominator(local_ele, ax, local_ind, d_N, d_OSEM);

#else // Implementation 3

				denominator_multi(local_ele, &axOSEM, read_imagef(d_OSEM, samplerIm, localInd).x);
				//denominator_multi(local_ele, &axOSEM, &d_OSEM[local_ind]);

#endif
				//if (d_N3 == 1)
				//	indX++;
				//else
				//	indY++;
				if (d_N3 == 1)
					localInd.x++;
				else
					localInd.y++;
				//local_ind += d_N3;
			}
#ifdef AF // Implementation 2

			nominator(MethodList, ax, local_sino, d_epsilon_mramla, d_epps, 1.f, d_sc_ra, idx);

#else // Implementation 3

			nominator_multi(&axOSEM, local_sino, d_epps, 1.f, d_sc_ra, idx);

#endif
			//local_ind = z_loop;
		}
#endif
#ifdef BP
		local_ind = tempk;
		for (uint ii = 0u; ii < d_N1; ii++) {
			//if (local_ind >= 1032192)
			//	continue;
#if defined(FP) && defined(BP)
			if (local_sino != 0.f) {
#endif
#ifdef AF
				rhs(MethodList, local_ele, ax, local_ind, d_N, d_output);
				//rhs(MethodList, local_ele, ax, local_ind, d_N, d_output);
#else
#ifdef ATOMIC // 64-bit atomics
				atom_add(&d_output[local_ind], convert_long(local_ele * axOSEM * TH));
				//atom_add(&d_output[local_ind], convert_long(local_ele * axOSEM * TH));
#elif defined(ATOMIC32)
				atomic_add(&d_output[local_ind], convert_int(local_ele * axOSEM * TH));
				//atomic_add(&d_output[local_ind], convert_int(local_ele* axOSEM* TH));
#else // 32-bit float atomics
				atomicAdd_g_f(&d_output[local_ind], (local_ele * axOSEM));
				//atomicAdd_g_f(&d_output[local_ind], (local_ele * axOSEM));
#endif
#endif
#if defined(FP) && defined(BP)
			}
#endif
			if (no_norm == 0u)
#ifdef ATOMIC // 64-bit atomics
				atom_add(&d_Summ[local_ind], convert_long(local_ele * TH));
			//atom_add(&d_Summ[local_ind], convert_long(local_ele * TH));
#elif defined(ATOMIC32)
				atomic_add(&d_Summ[local_ind], convert_int(local_ele * TH));
			//atomic_add(&d_Summ[local_ind], convert_int(local_ele* TH));
#else // 32-bit float atomics
				atomicAdd_g_f(&d_Summ[local_ind], local_ele);
			//atomicAdd_g_f(&d_Summ[local_ind], local_ele);
#endif
			local_ind += d_N3;
			//local_ind += d_N3;
		}
#endif
//		}
//		else {
//			//local_ele = templ_ijk;
//			//local_ind = z_loop;
//			for (uint ii = 0u; ii < d_N1; ii++) {
//#ifdef ATOMIC // 64-bit atomics
//				atom_add(&d_Summ[local_ind], convert_long(local_ele * TH));
//				//atom_add(&d_Summ[local_ind], convert_long(local_ele * TH));
//#elif defined(ATOMIC32)
//				atomic_add(&d_Summ[local_ind], convert_int(local_ele * TH));
//				//atomic_add(&d_Summ[local_ind], convert_int(local_ele* TH));
//#else // 32-bit float atomics
//				atomicAdd_g_f(&d_Summ[local_ind], local_ele);
//				//atomicAdd_g_f(&d_Summ[local_ind], local_ele);
//#endif
//				local_ind += d_N3;
//				//local_ind += d_N3;
//			}
//		}
#endif
#endif

#elif defined(ORTH) // Orthogonal or volume-based

#ifdef AF
orthDistancePerpendicularMulti3D(xcenter, center2, z_center, &temp, ax,
#elif defined(MBSREM)
orthDistancePerpendicularMulti3D(xcenter, center2, z_center, &temp, &axCOSEM,
#else
orthDistancePerpendicularMulti3D(xcenter, center2, z_center, &temp, &axOSEM,
#endif
	d_b, dd, d_db, d_N0, d_N1, d_Nxyz.z, tempk, local_norm, local_sino, d_N2, d_N3,
#if (defined(FP) && !defined(BP)) || defined(MBSREM)
	d_OSEM, 
#endif
	s, diff, kerroin, d_Nxy, no_norm, d_Summ, true, false, global_factor, bmin, bmax, Vmax, V,
#ifdef MBSREM
	MethodListOpenCL, d_alku, & axCOSEM, d_E, d_co, d_aco, & minimi, MBSREM_prepass, d_sc_ra, d_Amin, d_ACOSEM_lhs, idx
#else
	d_output, d_N, MethodList
#endif
#if !defined(CT) && defined(ATN)
	, d_atten
#endif
);

//#ifndef AF
#if defined(FP) && !defined(BP)
		//if (fp == 1) {
#ifdef AF // Implementation 2

			nominator(MethodList, ax, local_sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx);
			forwardProjectAF(d_output, ax, idx, d_N);

#else // Implementation 3
			nominator_multi(&axOSEM, local_sino, d_epps, temp, d_sc_ra, idx);
			d_output[idx] = axOSEM;
#endif
			return;
		//}
#endif
//#endif
#ifdef MBSREM
		if (d_alku == 0 && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1
			|| MethodListOpenCL.OSLCOSEM > 0) && local_sino != 0.f) {
			nominator_cosem(&axCOSEM, local_sino, d_epps, temp, d_sc_ra, idx);
		}

#else
#if defined(FP) && defined(BP) // Forward projection
		if (local_sino != 0.f) {

#ifdef AF
			nominator(MethodList, ax, local_sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx);
#else
			nominator_multi(&axOSEM, local_sino, d_epps, temp, d_sc_ra, idx);
#endif
		}
#endif
#endif
#ifdef BP
#ifdef AF
		orthDistancePerpendicularMulti3D(xcenter, center2, z_center, &temp, ax,
#elif defined(MBSREM)
		orthDistancePerpendicularMulti3D(xcenter, center2, z_center, &temp, &axCOSEM,
#else
		orthDistancePerpendicularMulti3D(xcenter, center2, z_center, &temp, &axOSEM,
#endif
			d_b, dd, d_db, d_N0, d_N1, d_Nxyz.z, tempk, local_norm, local_sino, d_N2, d_N3,
#if (defined(FP) && !defined(BP)) || defined(MBSREM)
			d_OSEM,
#endif
			s, diff, kerroin, d_Nxy, no_norm, d_Summ, false, true, global_factor, bmin, bmax, Vmax, V,
#ifdef MBSREM
			MethodListOpenCL, d_alku, & axCOSEM, d_E, d_co, d_aco, & minimi, MBSREM_prepass, d_sc_ra, d_Amin, d_ACOSEM_lhs, idx
#else
			d_output, d_N, MethodList
#endif
#if !defined(CT) && defined(ATN)
			, d_atten
#endif
		);

#endif
#endif
#endif

	}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///*
	else {
	//return;
		//float3 tu = { 0.f, 0.f, 0.f };
		//float3 t0 = { 1e8f, 1e8f, 1e8f };
		//float tc = 0.f;
		//float txu = 0.f, tyu = 0.f, tzu = 0.f, tc = 0.f, tx0 = 1e8f, ty0 = 1e8f, tz0 = 1e8f;
		float txu = 0.f, tyu = 0.f, tzu = 0.f, tc = 0.f, tx0 = 1e8f, ty0 = 1e8f, tz0 = 1e8f;
		bool skip = false, XY = true;

		// If the measurement is on a same ring
			// Z-coordinate (ring)
		if (fabs(diff.z) < 1e-6f) {
			//return;
			//tempk = convert_int((zs / d_zmax) * (d_NSlices - 1.f));
			tempk = convert_int(fabs(s.z - b.z) / d_d.z);
			skip = siddon_pre_loop_2D(b.x, b.y, diff.x, diff.y, d_max.x, d_max.y, d_d.x, d_d.y, d_Nxyz.x, d_Nxyz.y, &tempi, &tempj, &txu, &tyu, &Np, TYPE,
				s.y, s.x, d.y, d.x, &tc, &ux, &uy, &tx0, &ty0, &XY);
			//if (i.z == 0 && i.y == 100 && i.x == 150) {
			//	printf("tempi = %d\n", tempi);
			//	printf("tempj = %d\n", tempj);
			//	printf("tempk = %d\n", tempk);
			//	printf("Np = %d\n", Np);
			//	printf("tx0 = %f\n", tx0);
			//	printf("ty0 = %f\n", ty0);
			//	printf("tz0 = %f\n", tz0);
			//	printf("tc = %f\n", tc);
			//	//printf("s.x = %f\n", s.x);
			//	//printf("s.y = %f\n", s.y);
			//	//printf("s.z = %f\n", s.z);
			//	//printf("d.x = %f\n", d.x);
			//	//printf("d.y = %f\n", d.y);
			//	//printf("d.z = %f\n", d.z);
			//	//printf("b.x = %f\n", b.x);
			//	//printf("b.y = %f\n", b.y);
			//	//printf("b.z = %f\n", b.z);
			//}
			//skip = siddon_pre_loop_2D(d_bx, d_by, x_diff, y_diff, d_maxxx, d_maxyy, d_dx, d_dy, d_Nx, d_Ny, &tempi, &tempj, &txu, &tyu, &Np, TYPE,
			//	ys, xs, yd, xd, &tc, &ux, &uy, &tx0, &ty0);
		}
		else if (fabs(diff.y) < 1e-6f) {
			//return;
			tempj = perpendicular_start(b.y, d.y, d_d.y, d_Nxyz.y);
			skip = siddon_pre_loop_2D(b.x, b.z, diff.x, diff.z, d_max.x, d_max.z, d_d.x, d_d.z, d_Nxyz.x, d_Nxyz.z, &tempi, &tempk, &txu, &tzu, &Np, TYPE,
				s.z, s.x, d.z, d.x, &tc, &ux, &uz, &tx0, &tz0, &XY);
			XY = true;
#ifndef PRECOMPUTE
			if (d.y > d_max.y || d.y < b.y)
				skip = true;
#endif
		}
		else if (fabs(diff.x) < 1e-6f) {
			//return;
			tempi = perpendicular_start(b.x, d.x, d_d.x, d_Nxyz.x);
			skip = siddon_pre_loop_2D(b.y, b.z, diff.y, diff.z, d_max.y, d_max.z, d_d.y, d_d.z, d_Nxyz.y, d_Nxyz.z, &tempj, &tempk, &tyu, &tzu, &Np, TYPE,
				s.z, s.y, d.z, d.y, &tc, &uy, &uz, &ty0, &tz0, &XY);
			XY = false;
//#ifndef SIDDON
//			int apu_tempi = tempi;
//			float apu_txu = txu;
//			float apu_tx0 = tx0;
//			float apu_xdiff = diff.x;
//			float apu_xs = s.x;
//			int apu_iu = ux;
//			ux = uy;
//			uy = apu_iu;
//			tempi = tempj;
//			tempj = apu_tempi;
//			txu = tyu;
//			tyu = apu_txu;
//			tx0 = ty0;
//			ty0 = apu_tx0;
//			diff.x = diff.y;
//			diff.y = apu_xdiff;
//			s.x = s.y;
//			s.y = apu_xs;
//			d_N0 = d_Nxyz.y;
//			d_N1 = d_Nxyz.x;
//			d_N2 = d_Nxyz.y;
//			d_N3 = 1u;
//#ifdef ORTH // Orthogonal or volume-based
//			ycenter = x_center;
//			xcenter = y_center;
//#endif
//#endif
#ifndef PRECOMPUTE
			if (d.x > d_max.x || d.x < b.x)
				skip = true;
#endif
		}
		else {
			//return;
			skip = siddon_pre_loop_3D(b, diff, d_max, d_d, d_Nxyz, &tempi, &tempj, &tempk, &txu, &tyu, &tzu, &Np, TYPE, s, d, &tc, &ux, &uy, &uz, &tx0, &ty0, &tz0, &XY);
			//skip = siddon_pre_loop_3D(d_bx, d_by, d_bz, x_diff, y_diff, z_diff, d_maxxx, d_maxyy, d_bzb, d_dx, d_dy, d_dz, d_Nx, d_Ny, d_Nz, &tempi, &tempj, &tempk, &tyu, &txu, &tzu,
			//	&Np, TYPE, ys, xs, yd, xd, zs, zd, &tc, &ux, &uy, &uz, &tx0, &ty0, &tz0);
		}

		////if (i.x == 30 && i.y == 86 && i.z == 543) {
		//if (i.x == 100 && i.y == 100 && i.z == 0) {
		//	printf("tc = %f\n", tc);
		//	//printf("tyu = %f\n", tyu);
		//	printf("tempi = %d\n", tempi);
		//	printf("tempj = %d\n", tempj);
		//	printf("tempk = %d\n", tempk);
		//	//printf("idx = %u\n", idx);
		//	printf("tx0 = %f\n", tx0);
		//	printf("ty0 = %f\n", ty0);
		//	printf("tz0 = %f\n", tz0);			
		//	//return;
		//}
#ifndef PRECOMPUTE // No precomputation step performed
		if (skip)
			return;
#endif
#ifdef FIND_LORS // Precomputation phase
		for (uint ii = 0u; ii < Np; ii++) {
			temp_koko++;
			if (tz0 < ty0 && tz0 < tx0) {
				tempk += uz;
				tz0 += tzu;
			}
			else if (ty0 < tx0) {
				tempj += uy;
				ty0 += tyu;
			}
			else {
				tempi += ux;
				tx0 += txu;
			}
			if (tempi < 0 || tempi >= d_Nxyz.x || tempj < 0 || tempj >= d_Nxyz.y || tempk < 0 || tempk >= d_Nxyz.z)
				break;
		}
		d_lor[idx] = temp_koko;
		return;
#else // Not the precomputation phase
		float temp = 0.f;
#if defined(SIDDON) // Siddon
		const float tx0_a = tx0, ty0_a = ty0, tz0_a = tz0;
		const int tempi_a = tempi, tempj_a = tempj, tempk_a = tempk;
#endif
#if defined(SIDDON)// Siddon
		//L = native_sqrt(x_diff * x_diff + z_diff * z_diff + y_diff * y_diff);
		L = length(diff);
		const float tc_a = tc;
#ifdef TOF
		TOFDis(diff, tc, L, &D, &DD);
#endif
#if defined(DEC) && defined(TOF) // Save intermediate TOF results
		__private float store_elements[DEC * NBINS];
#elif defined(TOF)
		__private float store_elements[1];
#endif
#endif
#if defined(ATN) && !defined(SIDDON)
		//L = native_sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);
		L = length(diff);
#endif
#ifdef ORTH // Orthogonal or volume-based
		int tempi_b, u_b;
		uint d_Nb;
		uint local_ind = 0u;
		float t0_b, tu_b, diff_b, s_b, dd_b;
		uint3 d_NN = { d_Nxyz.x, d_Nxyz.y, d_Nxyz.z };
		//__local float* xCenter;
		//__local float* yCenter;
#ifdef DEC // Save intermediate results
		__private float store_elements[DEC];
		__private uint store_indices[DEC];
#else
		__private float store_elements[1];
		__private uint store_indices[1];
#endif
		//if (i.z == 543 && i.y == 4 && i.x == 115) {
		//	printf("tempi = %d\n", tempi);
		//	printf("tempj = %d\n", tempj);
		//	printf("tempk = %d\n", tempk);
		//	printf("Np = %d\n", Np);
		//		printf("tx0 = %f\n", tx0);
		//		printf("ty0 = %f\n", ty0);
		//		printf("tz0 = %f\n", tz0);		
		//		printf("tc = %f\n", tc);
		//		printf("xyvar = %d\n", xyvar);
		//			//printf("s.x = %f\n", s.x);
		//			//printf("s.y = %f\n", s.y);
		//			//printf("s.z = %f\n", s.z);
		//			//printf("d.x = %f\n", d.x);
		//			//printf("d.y = %f\n", d.y);
		//			//printf("d.z = %f\n", d.z);
		//			//printf("b.x = %f\n", b.x);
		//			//printf("b.y = %f\n", b.y);
		//			//printf("b.z = %f\n", b.z);
		//}
		//bool XY = true;
		if (!XY) {
		//if (fabs(diff.y) > fabs(diff.x)) {
		//if (fabs(s.x) > fabs(s.y)) {
		//if (ty0 > tx0) {
		//if ((tempj == d_Nxyz.y - 1 || tempj == 0) && tempi != 0 && tempi != d_Nxyz.x - 1) {
			tempi_b = tempi;
			tempi = tempj;
			tempj = tempi_b;
			s_b = s.x;
			s.x = s.y;
			s.y = s_b;
			dd_b = d.x;
			d.x = d.y;
			d.y = dd_b;
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
			xcenter = y_center;
			ycenter = x_center;
			d_NN.x = d_Nxyz.y;
			d_NN.y = d_Nxyz.x;
			d_N0 = d_Nxyz.y;
			d_N1 = d_Nxyz.x;
			d_N2 = d_Nxyz.x;
			d_N3 = 1;
			//XY = false;
		}
		else {
			//return;
			//if (tempj == d_Nxyz.y - 1 || tempj == 0) {
			//	tempi_b = tempi;
			//	tempi = tempj;
			//	tempj = tempi_b;
			//	s_b = s.x;
			//	s.x = s.y;
			//	s.y = s_b;
			//	dd_b = d.x;
			//	d.x = d.y;
			//	d.y = dd_b;
			//	u_b = ux;
			//	ux = uy;
			//	uy = u_b;
			//	t0_b = tx0;
			//	tx0 = ty0;
			//	ty0 = t0_b;
			//	tu_b = txu;
			//	txu = tyu;
			//	tyu = tu_b;
			//	diff_b = diff.x;
			//	diff.x = diff.y;
			//	diff.y = diff_b;
			//	xcenter = y_center;
			//	ycenter = x_center;
			//	d_NN.x = d_Nxyz.y;
			//	d_NN.y = d_Nxyz.x;
			//	d_N0 = d_Nxyz.y;
			//	d_N1 = d_Nxyz.x;
			//	d_N2 = d_Nxyz.x;
			//	d_N3 = 1;
			//	XY = false;
			//}
			//else {
				xcenter = x_center;
				ycenter = y_center;
				d_N0 = d_Nxyz.x;
				d_N1 = d_Nxyz.y;
				d_N2 = 1;
				d_N3 = d_Nxyz.x;
			//}
			//if (i.x == 80 && i.y == 6 && i.z == 20) {
			//	printf("tempi = %d\n", tempi);
			//	printf("tempj = %d\n", tempj);
			//	printf("tempk = %d\n", tempk);
			//}
		}
#if !defined(DEC)
		const float tx0_a = tx0, ty0_a = ty0, tz0_a = tz0;
		const int tempi_a = tempi, tempj_a = tempj, tempk_a = tempk;
#endif
		//if (tempi != 0 && tempi != d_Nxyz.x - 1 && i.z == 0 && i.y == 100 && i.x == 150) {
		//	printf("tempi = %d\n", tempi);
		//	printf("tempj = %d\n", tempj);
		//	printf("tempk = %d\n", tempk);
		//	//printf("i.x = %d\n", i.x);
		//	//printf("i.y = %d\n", i.y);
		//	printf("Np = %d\n", Np);
		//	printf("tx0 = %f\n", tx0);
		//	printf("ty0 = %f\n", ty0);
		//	printf("tz0 = %f\n", tz0);
		//	printf("tc = %f\n", tc);
		//	printf("s.x = %f\n", s.x);
		//	printf("s.y = %f\n", s.y);
		//	printf("s.z = %f\n", s.z);
		//	printf("d.x = %f\n", d.x);
		//	printf("d.y = %f\n", d.y);
		//	printf("d.z = %f\n", d.z);
		//	printf("diff.x = %f\n", diff.x);
		//	printf("diff.y = %f\n", diff.y);
		//}
		//bool ekaz = false;
		//bool ekay = false;
		//int temp1 = tempj;
		//int temp2 = tempk;
		for (uint ii = 0u; ii < Np; ii++) {
			if (tz0 < ty0 && tz0 < tx0) {
#ifdef ATN
				compute_attenuation(&tc, &jelppi, L, tz0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
#endif
				tempk += uz;
				tc = tz0;
				tz0 += tzu;
			}
			else if (ty0 < tx0) {
#ifdef ATN
				compute_attenuation(&tc, &jelppi, L, ty0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
#endif
				tempj += (uy);
				tc = tz0;
				ty0 += tyu;
			}
			else {
#ifdef ATN
				compute_attenuation(&tc, &jelppi, L, tx0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
#endif
				orthDistance3D(tempi, diff.y, diff.x, diff.z, xcenter[tempi], ycenter, z_center, &temp, tempj, tempk, local_sino,
#if (defined(FP) && !defined(BP)) || defined(MBSREM)
					d_OSEM,
#endif
					s.x, s.y, s.z, d_Nxy, kerroin, no_norm, RHS, SUMMA, d_Summ, d_N1, d_N2, d_N3, d_Nxyz.z, bmin, bmax, Vmax, V, 
					store_elements, XY, i, 
#ifdef AF
#ifdef MBSREM
					store_indices, &local_ind, &axACOSEM, MethodListOpenCL, d_E, idx, d_co, d_aco, &minimi, MBSREM_prepass, &axCOSEM, d_alku);
#else
					store_indices, & local_ind, d_output, d_N, MethodList, ax);
#endif
#else
					store_indices, & local_ind, d_output, d_N, MethodList, & axOSEM);
#endif
				//temp1 = tempj;
				//temp2 = tempk;
				tempi += ux;
				tc = tx0;
				tx0 += txu;
			}
#ifndef PRECOMPUTE
#ifndef DEC
			Np_n++;
#endif
			//if (i.x == 30 && i.y == 86 && i.z == 543) {
			//	printf("tempi = %d\n", tempi);
			//	printf("tempj = %d\n", tempj);
			//	printf("tempk = %d\n", tempk);
			//	printf("s.x = %f\n", s.x);
			//	printf("s.y = %f\n", s.y);
			//	printf("s.z = %f\n", s.z);
			//	printf("d.x = %f\n", d.x);
			//	printf("d.y = %f\n", d.y);
			//	printf("d.z = %f\n", d.z);
			//}
			if (tempi < 0 || tempi >= d_NN.x || tempj < 0 || tempj >= d_NN.y || tempk < 0 || tempk >= d_NN.z) {
				break;
			}
#endif
		}
		if (tempi >= 0 && tempi < d_NN.x) {
			int aloitus = tempi;
			//if (tempi < 0)
			//	aloitus = 0;
			//else if (aloitus >= d_NN.x)
			//	aloitus = d_NN.x - 1;
			if (tempj < 0)
				tempj = 0;
			else if (tempj >= d_NN.y)
				tempj = d_NN.y - 1;
			if (tempk < 0)
				tempk = 0;
			else if (tempk >= d_NN.z)
				tempk = d_NN.z - 1;
			if (ux < 0) {
				//const int minimi = max(aloitus - NSTEPS, 0);
				const int minimi = 0;
				for (int kk = aloitus; kk >= minimi; kk--) {
					orthDistance3D(kk, diff.y, diff.x, diff.z, xcenter[kk], ycenter, z_center, &temp, tempj, tempk, local_sino,
#if (defined(FP) && !defined(BP)) || defined(MBSREM)
						d_OSEM,
#endif
						s.x, s.y, s.z, d_Nxy, kerroin, no_norm, RHS, SUMMA, d_Summ, d_N1, d_N2, d_N3, d_Nxyz.z, bmin, bmax, Vmax, V,
						store_elements, XY, i,
#ifdef AF
#ifdef MBSREM
						store_indices, &local_ind, &axACOSEM, MethodListOpenCL, d_E, idx, d_co, d_aco, &minimi, MBSREM_prepass, &axCOSEM, d_alku);
#else
						store_indices, & local_ind, d_output, d_N, MethodList, ax);
#endif
#else
						store_indices, & local_ind, d_output, d_N, MethodList, & axOSEM);
#endif
				}
			}
			else {
				//const int maksimi = min(aloitus + NSTEPS, convert_int(d_N0));
				const int maksimi = convert_int(d_N0);
				for (int kk = aloitus; kk < maksimi; kk++) {
					orthDistance3D(kk, diff.y, diff.x, diff.z, xcenter[kk], ycenter, z_center, &temp, tempj, tempk, local_sino,
#if (defined(FP) && !defined(BP)) || defined(MBSREM)
						d_OSEM,
#endif
						s.x, s.y, s.z, d_Nxy, kerroin, no_norm, RHS, SUMMA, d_Summ, d_N1, d_N2, d_N3, d_Nxyz.z, bmin, bmax, Vmax, V,
						store_elements, XY, i,
#ifdef AF
#ifdef MBSREM
						store_indices, &local_ind, &axACOSEM, MethodListOpenCL, d_E, idx, d_co, d_aco, &minimi, MBSREM_prepass, &axCOSEM, d_alku);
#else
						store_indices, & local_ind, d_output, d_N, MethodList, ax);
#endif
#else
						store_indices, & local_ind, d_output, d_N, MethodList, & axOSEM);
#endif
				}
			}
		}
		//if (i.x == 126 && i.y == 35 && i.z == 697) {
		//	printf("axOSEM = %f\n", axOSEM);
		//}
//#ifdef PRECOMPUTE
//#ifdef CRYSTZ
//		if (xyz < 3 && fabs(z_diff) >= 1e-6f) {
//			if (xyz == 1)
//				tempi -= ux;
//			else if (xyz == 2)
//				tempj -= uy;
//			if ((tempk >= (d_Nz - 1) && uz > 0) || (tempk <= 0 && uz < 0) || uz == 0) {}
//			else {
//				tempk += uz;
//				alku = d_Nz;
//				loppu = 0;
//				if (uz > 0) {
//					loppu = tempk;
//				}
//				else if (uz < 0) {
//					alku = tempk + 1;
//				}
//				orthDistance3D(tempi, d_N0, d_N4, y_diff, x_diff, z_diff, y_center, x_center, z_center, &temp, d_N2, tempj, tempk,
//					local_sino, d_OSEM, xs, ys, zs, d_Nxy, kerroin, no_norm, RHS, SUMMA, d_Summ, d_N1, d_N3, alku,
//					ux, uy, loppu, bmin, bmax, Vmax, V, store_elements,
//#ifdef AF
//#ifdef MBSREM
//					store_indices, &local_ind, &axACOSEM, MethodListOpenCL, d_E, idx, d_co, d_aco, &minimi, MBSREM_prepass, &axCOSEM, d_alku);
//#else
//					store_indices, & local_ind, d_output, d_N, MethodList, ax);
//#endif
//#else
//					store_indices, & local_ind, d_output, d_N, MethodList, & axOSEM);
//#endif
//			}
//		}
//#endif
//#endif
#elif defined(SIDDON)
#if !defined(CT) || (defined(FP) && defined(CT))
		for (uint ii = 0U; ii < Np; ii++) {
			//local_ind = compute_ind(tempj, tempi, tempk, d_N0, d_N1, d_N, d_N3, d_Nxy);
			localInd = (int4)(tempi, tempj, tempk, 0);
			if (tz0 < ty0 && tz0 < tx0) {
				local_ele = compute_element(&tz0, &tc, L, tzu, uz, &tempk, &temp);
			}
			else if (ty0 < tx0) {
				local_ele = compute_element(&ty0, &tc, L, tyu, uy, &tempj, &temp);
			}
			else {
				local_ele = compute_element(&tx0, &tc, L, txu, ux, &tempi, &temp);
			}
#ifdef MATRIX
			local_ele /= L;
#ifdef ATN
			local_ele *= native_exp(local_ele * -read_imagef(d_atten, samplerIm, localInd).x));
#endif
#ifdef NORM
			local_ele *= local_norm;
#endif
#ifdef SCATTER
			local_ele *= local_scat;
#endif
			local_ele *= global_factor;
			indices[idx * maxLOR + ii] = local_ind;
			values[idx * maxLOR + ii] = local_ele;
#else
#ifdef ATN
			jelppi += (local_ele * -read_imagef(d_atten, samplerIm, localInd).x);
#endif
#if defined(TOF) && (defined(DEC) || defined(FP))
			const float TOFSum = TOFLoop(DD, local_ele, store_elements, TOFCenter, sigma_x, &D, ii * NBINS, d_epps);
#endif
#ifdef FP
#ifdef BP
			if (local_sino != 0.f) {
#endif
#ifdef AF // Implementation 2

#ifdef MBSREM
#ifdef TOF
				if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0u)
					denominatorTOF(ax, local_ele, d_OSEM, localInd.x, localInd.y, localInd.z, TOFSum, store_elements, DD, TOFCenter, sigma_x, &D, ii * NBINS, d_epps, d_N);
#else
				if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0u)
					denominator_multi(local_ele, &axCOSEM, read_imagef(d_OSEM, samplerIm, localInd).x);
#endif
#else
#ifdef TOF
				denominatorTOF(ax, local_ele, d_OSEM, localInd.x, localInd.y, localInd.z, TOFSum, store_elements, DD, TOFCenter, sigma_x, &D, ii * NBINS, d_epps, d_N);
#else
				denominator(local_ele, ax, localInd, d_N, d_OSEM);
#endif

#endif

#else // Implementation 3
#ifdef TOF
				denominatorTOF(ax, local_ele, d_OSEM, localInd, TOFSum, store_elements, DD, TOFCenter, sigma_x, &D, ii* NBINS, d_epps, d_N);
#else
				denominator_multi(local_ele, &axOSEM, read_imagef(d_OSEM, samplerIm, localInd).x);
#endif

#endif
#ifdef BP
			}
#endif
#endif
#endif
			if (i.x == 100 && i.y == 100 && i.z == 0) {
				printf("axOSEM = %f\n", axOSEM);
				printf("tempi = %d\n", tempi);
				printf("tempj = %d\n", tempj);
				printf("tempk = %d\n", tempk);
			}
#ifndef PRECOMPUTE
			Np_n++;
			if (tempi < 0 || tempi >= d_Nxyz.x || tempj < 0 || tempj >= d_Nxyz.y || tempk < 0 || tempk >= d_Nxyz.z) {
			//if (any(localInd < 0) || any(localInd > as_int4(d_Nxyz))) {
				break;
			}
#endif
		}
		//temp = L;
#else
		Np_n = Np;
#endif
#ifdef MATRIX
		rowInd[idx + 1] = Np_n;
		return;
#endif
#endif

#ifndef CT
		temp = 1.f / temp;
#ifdef ATN
		temp *= native_exp(jelppi);
#endif
#ifdef NORM
		temp *= local_norm;
#endif
#ifdef SCATTER
		temp *= d_scat[idx];
#endif
		temp *= global_factor;
#endif

#if defined(SIDDON) || !defined(DEC)
		tx0 = tx0_a, ty0 = ty0_a, tz0 = tz0_a;
		tempi = tempi_a, tempj = tempj_a, tempk = tempk_a;
#endif
#if defined(SIDDON)
		tc = tc_a;
#ifdef TOF
		D = DD;
#endif
//#else
//#if !defined(DEC)
//		xyz = 0u;
//#endif
#endif
#ifdef MBSREM
		if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && local_sino != 0.f && d_alku == 0u) {
#ifdef TOF
#pragma unroll NBINS
			for (int to = 0; to < NBINS; to++) {
				ax[to] *= temp;
				if (ax[to] < d_epps)
					ax[to] = d_epps;
#ifdef RANDOMS
				ax[to] += d_sc_ra[idx];
#endif
				ax[to] = d_Sino[idx + to * m_size] / ax[to];
			}
#else
#ifndef CT
			axCOSEM *= temp;
			if (axCOSEM < d_epps)
				axCOSEM = d_epps;
#endif
#ifdef RANDOMS
			axCOSEM += d_sc_ra[idx];
#endif
#ifdef CT
			axCOSEM = native_exp(-axCOSEM) / local_sino;
#else
			axCOSEM = local_sino / axCOSEM;
#endif
#endif
		}
		RHS = true;
#else
//#ifndef AF
#if defined(FP) && !defined(BP)
		//if (fp == 1) {
#ifdef TOF
			nominatorTOF(MethodList, ax, d_Sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx, m_size, local_sino);
#pragma unroll NBINS
			for (int to = 0; to < NBINS; to++)
				d_output[idx + to * m_size] = ax[to];
#else
			//if (i.x == 89 && i.y == 89 && i.z == 20) {
			//	printf("idx = %d\n", idx);
			//	printf("axOSEM = %f\n", axOSEM);
			//}
			//if (i.x == 82 && i.y == 82 && i.z == 20) {
			//	printf("idx = %d\n", idx);
			//	printf("axOSEM = %f\n", axOSEM);
			//}
#ifdef AF // Implementation 2

			nominator(MethodList, ax, local_sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx);
			forwardProjectAF(d_output, ax, idx, d_N);

#else // Implementation 3
			nominator_multi(&axOSEM, local_sino, d_epps, temp, d_sc_ra, idx);
			d_output[idx] = axOSEM;
#endif
#endif
			return;
		//}
#endif
//#endif

#if defined(FP) && defined(BP)
		if (local_sino != 0.f) {
#ifdef TOF
			nominatorTOF(MethodList, ax, d_Sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx, m_size, local_sino);
#else
#ifdef AF
			nominator(MethodList, ax, local_sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx);
#else
			nominator_multi(&axOSEM, local_sino, d_epps, temp, d_sc_ra, idx);
#endif
#endif
		}
#endif
			RHS = true;
		//else
		//	SUMMA = true;
#endif
#if defined(BP) || defined(MBSREM)
		// Add additional computations before backprojection here
#ifdef ORTH
#ifdef DEC
		for (uint ii = 0u; ii < local_ind; ii++) {
#ifdef CT
			const float local_ele = store_elements[ii];
#else
			const float local_ele = store_elements[ii] * temp;
#endif
			const uint localind = store_indices[ii];
			if (RHS) {
#ifdef AF
				rhs(MethodList, local_ele, ax, localind, d_N, d_output);
#else
#ifdef ATOMIC
				atom_add(&d_output[localind], convert_long(local_ele * axOSEM * TH));
#elif defined(ATOMIC32)
				atomic_add(&d_output[localind], convert_int(local_ele * axOSEM * TH));
#else
				atomicAdd_g_f(&d_output[localind], (local_ele * axOSEM));
#endif
#endif
			}
			if (no_norm == 0u)
#ifdef ATOMIC
				atom_add(&d_Summ[localind], convert_long(local_ele* TH));
#elif defined(ATOMIC32)
				atomic_add(&d_Summ[localind], convert_int(local_ele * TH));
#else
				atomicAdd_g_f(&d_Summ[localind], local_ele);
#endif
		}
#else
			//if (!XY) {
			//	tempi_b = tempi;
			//	tempi = tempj;
			//	tempj = tempi_b;
			//	t0_b = tx0;
			//	tx0 = ty0;
			//	ty0 = t0_b;
			//}
		for (uint ii = 0u; ii < Np_n; ii++) {
			if (tz0 < ty0 && tz0 < tx0) {
				tempk += uz;
				tz0 += tzu;
			}
			else if (ty0 < tx0) {
				tempj += (uy);
				ty0 += tyu;
			}
			else {
				orthDistance3D(tempi, diff.y, diff.x, diff.z, xcenter[tempi], ycenter, z_center, &temp, tempj, tempk, local_sino,
#if (defined(FP) && !defined(BP)) || defined(MBSREM)
					d_OSEM,
#endif
					s.x, s.y, s.z, d_Nxy, kerroin, no_norm, RHS, SUMMA, d_Summ, d_N1, d_N2, d_N3, d_Nxyz.z, bmin, bmax, Vmax, V,
					store_elements, XY, i,
#ifdef AF
#ifdef MBSREM
					store_indices, &local_ind, &axACOSEM, MethodListOpenCL, d_E, idx, d_co, d_aco, &minimi, MBSREM_prepass, &axCOSEM, d_alku);
#else
					store_indices, & local_ind, d_output, d_N, MethodList, ax);
#endif
#else
					store_indices, & local_ind, d_output, d_N, MethodList, & axOSEM);
#endif
				tempi += ux;
				tx0 += txu;
			}
		}
		if (tempi >= 0 && tempi < d_NN.x) {
			int aloitus = tempi;
			//if (tempi < 0)
			//	aloitus = 0;
			//else if (aloitus >= d_NN.x)
			//	aloitus = d_NN.x - 1;
			if (tempj < 0)
				tempj = 0;
			else if (tempj >= d_NN.y)
				tempj = d_NN.y - 1;
			if (tempk < 0)
				tempk = 0;
			else if (tempk >= d_NN.z)
				tempk = d_NN.z - 1;
			if (ux < 0) {
				//const int minimi = max(tempi - NSTEPS, 0);
				const int minimi = 0;
				for (int kk = tempi - 1; kk >= minimi; kk--) {
					orthDistance3D(kk, diff.y, diff.x, diff.z, xcenter[kk], ycenter, z_center, &temp, tempj, tempk, local_sino,
#if (defined(FP) && !defined(BP)) || defined(MBSREM)
						d_OSEM,
#endif
						s.x, s.y, s.z, d_Nxy, kerroin, no_norm, RHS, SUMMA, d_Summ, d_N1, d_N2, d_N3, d_Nxyz.z, bmin, bmax, Vmax, V,
						store_elements, XY, i,
#ifdef AF
#ifdef MBSREM
						store_indices, &local_ind, &axACOSEM, MethodListOpenCL, d_E, idx, d_co, d_aco, &minimi, MBSREM_prepass, &axCOSEM, d_alku);
#else
						store_indices, & local_ind, d_output, d_N, MethodList, ax);
#endif
#else
						store_indices, & local_ind, d_output, d_N, MethodList, & axOSEM);
#endif
				}
			}
			else {
				//const int maksimi = min(tempi + NSTEPS, d_N0);
				const int maksimi = convert_int(d_N0);
				for (int kk = tempi + 1; kk < maksimi; kk++) {
					orthDistance3D(kk, diff.y, diff.x, diff.z, xcenter[kk], ycenter, z_center, &temp, tempj, tempk, local_sino,
#if (defined(FP) && !defined(BP)) || defined(MBSREM)
						d_OSEM,
#endif
						s.x, s.y, s.z, d_Nxy, kerroin, no_norm, RHS, SUMMA, d_Summ, d_N1, d_N2, d_N3, d_Nxyz.z, bmin, bmax, Vmax, V,
						store_elements, XY, i,
#ifdef AF
#ifdef MBSREM
						store_indices, &local_ind, &axACOSEM, MethodListOpenCL, d_E, idx, d_co, d_aco, &minimi, MBSREM_prepass, &axCOSEM, d_alku);
#else
						store_indices, & local_ind, d_output, d_N, MethodList, ax);
#endif
#else
						store_indices, & local_ind, d_output, d_N, MethodList, & axOSEM);
#endif
				}
			}
		}
#ifdef MBSREM
		if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1 && d_alku == 0u)
			d_Amin[idx] = minimi;
		if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u) {
#ifdef RANDOMS
				axACOSEM += d_sc_ra[idx];
#endif
			d_ACOSEM_lhs[idx] = axACOSEM;
		}
#endif

#endif
#elif defined(SIDDON)
		//if (RHS) {
			for (uint ii = 0u; ii < Np_n; ii++) {
				local_ind = compute_ind(tempj, tempi, tempk, d_N0, d_N1, d_N, d_N3, d_Nxy);
#if defined(TOF) || defined(MBSREM)
				localInd = (int4)(tempi, tempj, tempk, 0);
#endif
				if (tz0 < ty0 && tz0 < tx0) {
					local_ele = compute_element_2nd(&tz0, &tc, L, tzu, uz, &tempk, temp);
				}
				else if (ty0 < tx0) {
					local_ele = compute_element_2nd(&ty0, &tc, L, tyu, uy, &tempj, temp);
				}
				else {
					local_ele = compute_element_2nd(&tx0, &tc, L, txu, ux, &tempi, temp);
				}
#ifdef TOF
#ifndef DEC
				const float TOFSum = TOFLoop(DD, local_ele / temp, store_elements, TOFCenter, sigma_x, &D, ii * NBINS, d_epps);
#endif
				backprojectTOF(local_ind, localInd, local_ele, ii * NBINS, store_elements, ax, d_Summ, local_sino,
#ifndef DEC
					temp, sigma_x, &D, DD, TOFCenter, d_epps, TOFSum,
#endif
#ifdef MBSREM
					MethodListOpenCL, d_alku, MBSREM_prepass, minimi, axACOSEM, d_OSEM, d_E, d_co, d_aco, idx, m_size);
#else
					d_output, no_norm, d_N);
#endif
#else
#ifdef MBSREM
				if (d_alku == 0u) {
					if (MBSREM_prepass == 1)
#ifdef ATOMIC
						atom_add(&d_Summ[local_ind], convert_long(local_ele * TH));
#elif defined(ATOMIC32)
						atomic_add(&d_Summ[local_ind], convert_int(local_ele* TH));
#else
						atomicAdd_g_f(&d_Summ[local_ind], local_ele);
#endif
					if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
						if (local_ele < minimi && local_ele > 0.f)
							minimi = local_ele;
						d_E[idx] += local_ele;
					}
					if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino != 0.f)
#ifdef ATOMIC
						atom_add(&d_co[local_ind], convert_long(axCOSEM * local_ele * TH));
#elif defined(ATOMIC32)
						atomic_add(&d_co[local_ind], convert_int(axCOSEM* local_ele* TH));
#else
						atomicAdd_g_f(&d_co[local_ind], axCOSEM * local_ele);
#endif
					if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino != 0.f)
#ifdef ATOMIC
						atom_add(&d_aco[local_ind], convert_long(axCOSEM * local_ele * TH));
#elif defined(ATOMIC32)
						atomic_add(&d_aco[local_ind], convert_int(axCOSEM* local_ele* TH));
#else
						atomicAdd_g_f(&d_aco[local_ind], axCOSEM * local_ele);
#endif
				}
				if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u)
					axACOSEM += (local_ele * read_imagef(d_OSEM, samplerIm, localInd).x);
					//axACOSEM += (local_ele * d_OSEM[local_ind]);
#else
				if (no_norm == 0u)
#ifdef ATOMIC
					atom_add(&d_Summ[local_ind], convert_long(local_ele * TH));
#elif defined(ATOMIC32)
					atomic_add(&d_Summ[local_ind], convert_int(local_ele* TH));
#else
					atomicAdd_g_f(&d_Summ[local_ind], local_ele);
#endif
#if defined(FP) && defined(BP)
				if (local_sino != 0.f) {
#endif
#ifdef AF
					rhs(MethodList, local_ele, ax, local_ind, d_N, d_output);
#else

#ifdef ATOMIC
					atom_add(&d_output[local_ind], convert_long(local_ele * axOSEM * TH));
#elif defined(ATOMIC32)
					atomic_add(&d_output[local_ind], convert_int(local_ele * axOSEM * TH));
#else
					atomicAdd_g_f(&d_output[local_ind], (local_ele * axOSEM));
#endif
#endif
#if defined(FP) && defined(BP)
				}
#endif
#endif
#endif
#if defined(CT) && !defined(FP)
				if (tempi < 0 || tempi >= d_Nxyz.x || tempj < 0 || tempj >= d_Nxyz.y || tempk < 0 || tempk >= d_Nxyz.z) {
					break;
				}
#endif
			}
#endif
		//}
//		else {
//			for (uint ii = 0u; ii < Np_n; ii++) {
//				local_ind = compute_ind(tempj, tempi, tempk, d_N0, d_N1, d_N, d_N3, d_Nxy);
//				if (tz0 < ty0 && tz0 < tx0) {
//					local_ele = compute_element_2nd(&tz0, &tc, L, tzu, uz, &tempk, temp);
//				}
//				else if (ty0 < tx0) {
//					local_ele = compute_element_2nd(&ty0, &tc, L, tyu, uy, &tempj, temp);
//				}
//				else {
//					local_ele = compute_element_2nd(&tx0, &tc, L, txu, ux, &tempi, temp);
//				}
//#ifdef TOF
//#ifndef DEC
//				const float TOFSum = TOFLoop(DD, local_ele / temp, store_elements, TOFCenter, sigma_x, &D, ii * NBINS, d_epps);
//#endif
//				sensTOF(local_ind, local_ele, ii* NBINS, store_elements, d_Summ, 
//#ifndef DEC
//					temp, sigma_x, & D, DD, TOFCenter, d_epps, TOFSum,
//#endif
//#ifdef MBSREM
//					MethodListOpenCL, d_alku, MBSREM_prepass, minimi, axACOSEM, d_OSEM, d_E, idx, m_size,
//#endif
//					no_norm);
//#else
//#ifdef MBSREM
//				if (d_alku == 0u) {
//					if (MBSREM_prepass == 1)
//#ifdef ATOMIC
//						atom_add(&d_Summ[local_ind], convert_long(local_ele * TH));
//#elif defined(ATOMIC32)
//						atomic_add(&d_Summ[local_ind], convert_int(local_ele* TH));
//#else
//						atomicAdd_g_f(&d_Summ[local_ind], local_ele);
//#endif
//					if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
//						if (local_ele < minimi && local_ele > 0.f)
//							minimi = local_ele;
//						d_E[idx] += local_ele;
//					}
//				}
//				if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u)
//					axACOSEM += (local_ele * d_OSEM[local_ind]);
//#else
//#ifdef ATOMIC
//				atom_add(&d_Summ[local_ind], convert_long(local_ele* TH));
//#elif defined(ATOMIC32)
//				atomic_add(&d_Summ[local_ind], convert_int(local_ele* TH));
//#else
//				atomicAdd_g_f(&d_Summ[local_ind], local_ele);
//#endif
//#endif
//#endif
//#if defined(CT) && !defined(FP)
//				if (tempj < 0 || tempi < 0 || tempk < 0 || tempi >= d_N0 || tempj >= d_N1 || tempk >= d_Nz) {
//					break;
//				}
//#endif
//			}
//		}
#ifdef MBSREM
#ifdef TOF
#pragma unroll NBINS
		for (long to = 0L; to < NBINS; to++) {
			if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1 && d_alku == 0u)
				d_Amin[idx + to * m_size] = minimi[to];
			if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u) {
#ifdef RANDOMS
				axACOSEM[to] += d_sc_ra[idx];
#endif
				d_ACOSEM_lhs[idx + to * m_size] = axACOSEM[to];
			}
		}
#else
		if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1 && d_alku == 0u)
			d_Amin[idx] = minimi;
		if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u) {
#ifdef RANDOMS
			axACOSEM += d_sc_ra[idx];
#endif
#ifdef CT
			d_ACOSEM_lhs[idx] = native_exp(-axACOSEM);
#else
			d_ACOSEM_lhs[idx] = axACOSEM;
#endif
		}
#endif
#endif
#endif
#endif
	}
//*/
}

