
/*******************************************************************************************************************************************
* A matrix free improved Siddon's for multiple rays. This function calculates Summ = sum(A,1) (sum of every row) and rhs = A*(y./(A'*x)), 
* where A is the system matrix, y the measurements and x the estimate/image.
*
* Used by implementations 2 and 3.
*
* This version goes through all the LORs and determines on-the-fly if they intersect with the voxel space. Uses (optionally) multiple rays.
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
* d_size_x = the number of detector elements,
* d_det_per_ring = number of detectors per ring,
* d_Nxy = d_Nx * d_Ny,
* fp = if 1, then only forward projection is computed, if 2 only backprojection, if 0 then both
* d_atten = attenuation data (images),
* d_norm = normalization coefficients,
* d_Summ = buffer for d_Summ,
* d_lor = number of pixels that each LOR traverses,
* d_pseudos = location of pseudo rings,
* d_x/y/z_det = detector x/y/z-coordinates,
* d_xy/zindex = for sinogram format they determine the detector
* indices corresponding to each sinogram bin (unused with raw data),
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
* d_output = rhs values for OSEM/MLEM,
* d_OSEM = OSEM/MLEM estimate,
* d_Summ = Sensitivity image
*
* [1] Jacobs, F., Sundermann, E., De Sutter, B., Christiaens, M. Lemahieu, I. (1998). A Fast Algorithm to Calculate the Exact Radiological 
* Path through a Pixel or Voxel Space. Journal of computing and information technology, 6 (1), 89-94.
*
* Copyright (C) 2019-2022  Ville-Veikko Wettenhovi
*
* This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it wiL be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************************************************************************/

// Matrix free Improved Siddon's algorithm
#if (defined(CT) || defined(PET) || defined(SPECT)) && !defined(LISTMODE)
__kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, 1)))
#else
__kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(LOCAL_SIZE, 1, 1)))
#endif
void proj1SiddonMultiRay(const float global_factor, const float d_epps, const uint d_N, const uint3 d_Nxyz, const float3 d_d,
	const float3 b, const float3 d_max, const uint d_size_x, const uint d_det_per_ring, const uint d_Nxy, 
	const float sigma_x, const float2 dc_z, 
#ifdef MASKFP
	__read_only image2d_t maskFP,
#endif
#ifdef TOF
	__constant float* TOFCenter,
#endif
#if !defined(CT) && defined(ATN)
	__read_only image3d_t d_atten,
#endif
	__constant uchar* MethodList, const float d_epsilon_mramla,
#if defined(CT) || defined(SPECT) || defined(PET)
	const uint d_sizey, const long d_nProjections,
#endif
#if defined(LISTMODE)
	const __global float* restrict d_xy, const __global float* restrict d_z,
#else
	__constant float* d_xy, __constant float* d_z,
#endif
	const __global float* restrict d_norm, const __global float* restrict d_scat, __global CAST* restrict d_Summ, const uchar fp,
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
#ifndef MBSREM
__global CAST* d_output, const uchar no_norm, const ulong m_size,
#else
const uint d_alku, const uchar MBSREM_prepass, __global float* d_ACOSEM_lhs, __global float* d_Amin, __global CAST* d_co,
__global CAST* d_aco, __global float* d_E, const ulong m_size, const RecMethodsOpenCL MethodListOpenCL,
#endif
const ulong cumsum) {
	// Get the current global index
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
#ifdef TOF
	float local_sino = 0.f;
#if defined(BP) && defined(FP)
#pragma unroll NBINS
	for (long to = 0L; to < NBINS; to++)
		local_sino += d_Sino[idx + m_size * to];
#endif
#else
#if defined(BP) && defined(FP)
	const float local_sino = (d_Sino[idx]);
#else
	const float local_sino = 0.f;
#endif
#endif
#if !defined(MBSREM) && defined(BP) && defined(FP)
	if (no_norm == 1u && local_sino == 0.f)
		return;
#elif defined(MBSREM)
	const uchar no_norm = 0u;
#endif

	bool RHS = true;
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
	//bool RHS = local_sino != 0.f ? true : false;
	float ax[NROLLS];
#pragma unroll
	for (uint kk = 0; kk < NROLLS; kk++)
		ax[kk] = 0.f;
#endif
#else
//#if defined(FP) && defined(BP)
//	bool RHS = local_sino != 0.f ? true : false;
//#else
//	bool RHS = true;
//#endif
#ifdef TOF
	float ax[NBINS];
#pragma unroll NBINS
	for (uint to = 0; to < NBINS; to++)
		ax[to] = 0.f;
#else
	float axOSEM = 0.f;
#endif
#endif
#ifdef BP
#ifndef AF
	if (fp == 2) {
#ifdef TOF
#pragma unroll NBINS
		for (uint to = 0; to < NBINS; to++)
			ax[to] = d_OSEM[idx + to * m_size + cumsum];
#else
		axOSEM = d_OSEM[idx + cumsum];
#endif
	}
#endif
#endif
#ifdef TOF
	float D = 0.f;
#endif
	uint d_N0 = d_Nxyz.x;
	uint d_N1 = d_Nxyz.y;
	uint d_N2 = 1u;
	uint d_N3 = d_Nxyz.x;
	float jelppi = 0.f;
	float temp = 0.f;
	int tempi_a[N_RAYS];
	int tempj_a[N_RAYS];
	int tempk_a[N_RAYS];
	uint d_N0_a[N_RAYS];
	int iu_a[N_RAYS];
	int ju_a[N_RAYS];
	int ku_a[N_RAYS];
	float tx0_a[N_RAYS];
	float ty0_a[N_RAYS];
	float tz0_a[N_RAYS];
	float tc_a[N_RAYS];
	float txu_a[N_RAYS];
	float tyu_a[N_RAYS];
	float tzu_a[N_RAYS];
	float LL[N_RAYS];
	uint Np_n[N_RAYS];
	bool pass[N_RAYS];
#ifdef TOF
	float DD[N_RAYS];
#if defined(DEC) // Save intermediate TOF results
	__private float store_elements[DEC * NBINS];
#else
	__private float store_elements[1];
#endif
#endif
	int lor = -1;
	// Load the next detector index
	// raw list-mode data
//#pragma unroll N_RAYS3D
	for (int lorZ = 0u; lorZ < N_RAYS3D; lorZ++) {
//#pragma unroll N_RAYS2D
		for (int lorXY = 0u; lorXY < N_RAYS2D; lorXY++) {
			lor++;
			float3 s, d;
#if defined(LISTMODE)
			getDetectorCoordinatesListmode(d_xy, &s, &d, idx, lorXY, lorZ, dc_z);
#elif defined(RAW) // raw data
			getDetectorCoordinatesRaw(d_xy, d_z, d_L, d_det_per_ring, idx, &s, &d, lorXY, lorZ, dc_z); // Sinogram data
#elif !defined(SUBSETS) // Precomputation phase
			getDetectorCoordinatesFullSinogram(d_size_x, i, &s, &d, d_xy, d_z, lorXY, lorZ, dc_z);
#else // Not the precomputation phase
			getDetectorCoordinates(d_xyindex, d_zindex, idx, &s, &d, d_xy, d_z, lorXY, lorZ, dc_z);
#endif
			//if (i.x == 100 && i.y == 100 && i.z == 0) {
			//	//printf("i.x = %d\n", i.x);
			//	//printf("i.y = %d\n", i.y);
			//	//printf("i.z = %d\n", i.z);
			//	//printf("d_size_x = %d\n", d_size_x);
			//	//printf("d_sizey = %d\n", d_sizey);
			//	printf("lorXY = %d\n", lorXY);
			//	printf("lorZ = %d\n", lorZ);
			//	printf("lor = %d\n", lor);
			//	//printf("N_RAYS2D = %d\n", N_RAYS2D);
			//	//printf("N_RAYS3D = %d\n", N_RAYS3D);
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
			//	//printf("s.x = % f\n", s.x);
			//	//printf("s.y = % f\n", s.y);
			//	//printf("s.z = % f\n", s.z);
			//	//printf("d.x = % f\n", d.x);
			//	//printf("d.y = % f\n", d.y);
			//	//printf("d.z = % f\n", d.z);
			//	//printf("axOSEM = % f\n", axOSEM);
			//}
			// Calculate the x, y and z distances of the detector pair
			float3 diff = d - s;
			if (all(diff == 0.f))
				return;
			pass[lor] = false;
			uint Np = 0u;
			Np_n[lor] = 0u;
			int4 localInd = { 0, 0, 0, 0 };
			// If the measurement is on a same ring
			if (fabs(diff.z) < 1e-6f && (fabs(diff.y) < 1e-6f || fabs(diff.x) < 1e-6f)) {
				float d_b, dd, d_db, d_d2;
				if (fabs(diff.y) < 1e-6f && d.y <= d_max.y && d.y >= b.y && s.y <= d_max.y && s.y >= b.y) {
					d_b = b.y;
					dd = d.y;
					d_db = d_d.y;
					d_d2 = d_d.x;
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
					tx0_a[lor] = 1e7f, ty0_a[lor] = 1e9f;
					pass[lor] = true;
				}
				else if (fabs(diff.x) < 1e-6f && d.x <= d_max.x && d.x >= b.x && s.x <= d_max.x && s.x >= b.x) {
					d_b = b.x;
					dd = d.x;
					d_d2 = d_d.y;
					d_db = d_d.x;
					tx0_a[lor] = 1e9f, ty0_a[lor] = 1e7f;
					pass[lor] = true;
				}
				else
					return;
				if (pass[lor]) {
					Np_n[lor] = d_N1;
					int tempk = convert_int(fabs(s.z - b.z) / d_d.z);
					localInd.z = tempk;
					//uint apu = 0u;
					float element = 0.f;
					perpendicular_elements(d_b, d_db, d_N0, dd, d_d2, d_N1, &element, &localInd, &tempk, d_N2, d_N3, d_norm, idx, global_factor,
						d_scat
#if !defined(CT) && defined(ATN)
						, d_atten
#endif
					);
					//const float element = perpendicular_elements_multiray(d_b, d_d, d_N0, dd, d_d2, d_N1, d_atten, &apu, tempk, d_N2, d_N3, &jelppi);
					temp += element;
					tempk_a[lor] = tempk;
#ifdef FP
#ifdef TOF
					float dI = (d_d2 * d_N1) / 2.f * -sign(y_diff);
					D = dI;
					DD[lor] = D;
					uint local_ind = apu;
					for (uint ii = 0; ii < d_N1; ii++) {
						const float TOFSum = TOFLoop(DD[lor], d_d2, store_elements, TOFCenter, sigma_x, &D, ii * NBINS, d_epps);
						//denominatorTOF(ax, d_d2, d_OSEM, local_ind, TOFSum, store_elements, DD[lor], TOFCenter, sigma_x, &D, ii * NBINS, d_epps, d_N);
						denominatorTOF(ax, d_d2, d_OSEM, localInd, TOFSum, store_elements, DD[lor], TOFCenter, sigma_x, &D, ii * NBINS, d_epps, d_N);
						if (d_N3 == 1)
							localInd.x++;
						else
							localInd.y++;
						//local_ind += d_N3;
					}
#else
#ifdef MBSREM
					if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && local_sino != 0.f && d_alku == 0u) {
						for (uint k = 0u; k < d_N1; k++)
							//axCOSEM += (d_d2 * d_OSEM[apu + k * d_N3]);
							axCOSEM += (d_d2 * read_imagef(d_OSEM, samplerIm, localInd).x);
						if (d_N3 == 1)
							localInd.x++;
						else
							localInd.y++;
					}
#else
					if (RHS) {
						for (uint k = 0u; k < d_N1; k++) {
#ifdef AF
							//denominator(d_d2, ax, apu + k * d_N3, d_N, d_OSEM);
							denominator(d_d2, ax, localInd, d_N, d_OSEM);
#else
							axOSEM += (d_d2 * read_imagef(d_OSEM, samplerIm, localInd).x);
#endif
							if (d_N3 == 1)
								localInd.x++;
							else
								localInd.y++;
						}
					}
#endif
#endif
#endif
				}
			}
			else {
				//pass[lor] = false;
				//return;
				int tempi = 0, tempj = 0, tempk = 0, ux= 0, uy = 0, uz = 0;
				float txu = 0.f, tyu = 0.f, tzu = 0.f, tc = 0.f, tx0 = 1e8f, ty0 = 1e8f, tz0 = 1e8f;
				bool skip = false, XY = true;
				if (fabs(diff.z) < 1e-6f) {
					tempk = convert_int(fabs(s.z - b.z) / d_d.z);
					skip = siddon_pre_loop_2D(b.x, b.y, diff.x, diff.y, d_max.x, d_max.y, d_d.x, d_d.y, d_Nxyz.x, d_Nxyz.y, &tempi, &tempj, &txu, &tyu, &Np, TYPE,
						s.y, s.x, d.y, d.x, &tc, &ux, &uy, &tx0, &ty0, &XY);
				}
				else if (fabs(diff.y) < 1e-6f) {
					tempj = perpendicular_start(b.y, d.y, d_d.y, d_Nxyz.y);
					skip = siddon_pre_loop_2D(b.x, b.z, diff.x, diff.z, d_max.x, d_max.z, d_d.x, d_d.z, d_Nxyz.x, d_Nxyz.z, &tempi, &tempk, &txu, &tzu, &Np, TYPE,
						s.z, s.x, d.z, d.x, &tc, &ux, &uz, &tx0, &tz0, &XY);
					if (d.y > d_max.y || d.y < b.y)
						skip = true;
				}
				else if (fabs(diff.x) < 1e-6f) {
					tempi = perpendicular_start(b.x, d.x, d_d.x, d_Nxyz.x);
					skip = siddon_pre_loop_2D(b.y, b.z, diff.y, diff.z, d_max.y, d_max.z, d_d.y, d_d.z, d_Nxyz.y, d_Nxyz.z, &tempj, &tempk, &tyu, &tzu, &Np, TYPE,
						s.z, s.y, d.z, d.y, &tc, &uy, &uz, &ty0, &tz0, &XY);
					if (d.x > d_max.x || d.x < b.x)
						skip = true;
				}
				else {
					skip = siddon_pre_loop_3D(b, diff, d_max, d_d, d_Nxyz, &tempi, &tempj, &tempk, &txu, &tyu, &tzu, &Np, TYPE, s, d, &tc, &ux, &uy, &uz, &tx0, &ty0, &tz0, &XY);
				}
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
				if (!skip) {
					pass[lor] = true;
					//return;
				}
				if (pass[lor]) {
					float apuvar = 0.f;
					LL[lor] = length(diff);
					//LL[lor] = native_sqrt(x_diff * x_diff + z_diff * z_diff + y_diff * y_diff);
					tx0_a[lor] = tx0, ty0_a[lor] = ty0, tz0_a[lor] = tz0, tc_a[lor] = tc;
					txu_a[lor] = txu, tyu_a[lor] = tyu, tzu_a[lor] = tzu;
					tempi_a[lor] = tempi, tempj_a[lor] = tempj, tempk_a[lor] = tempk;
					d_N0_a[lor] = d_N0;
					//uint temp_ijk = compute_ind(tempj, tempi, tempk, d_N0, d_N1, d_N, d_N3, d_Nxy);
					//temp_ijk_a[lor] = temp_ijk;
					iu_a[lor] = ux, ju_a[lor] = uy, ku_a[lor] = uz;
#ifdef TOF
					TOFDis(diff, tc, L, &D, &DD);
#endif
					float local_ele;
					//if (i.x == 100 && i.y == 100 && i.z == 0) {
					//	printf("tempi = %d\n", tempi);
					//	printf("tempj = %d\n", tempj);
					//	printf("tempk = %d\n", tempk);
					//}
					for (uint ii = 0u; ii < Np; ii++) {
						//const uint local_ind = compute_ind(tempj, tempi, tempk, d_N0, d_N1, d_N, d_N0, d_Nxy);
						localInd = (int4)(tempi, tempj, tempk, 0);
						if (tz0 < ty0 && tz0 < tx0) {
							local_ele = compute_element(&tz0, &tc, LL[lor], tzu, uz, &tempk, &temp);
						}
						else if (ty0 < tx0) {
							local_ele = compute_element(&ty0, &tc, LL[lor], tyu, uy, &tempj, &temp);
						}
						else {
							local_ele = compute_element(&tx0, &tc, LL[lor], txu, ux, &tempi, &temp);
						}
#ifdef ATN
						jelppi += (local_ele * -read_imagef(d_atten, samplerIm, localInd).x);
#endif
#ifdef TOF
						const float TOFSum = TOFLoop(DD[lor], local_ele, store_elements, TOFCenter, sigma_x, &D, ii * NBINS, d_epps);
#endif
#ifdef FP
#ifdef MBSREM
#ifdef TOF
						if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && local_sino != 0.f && d_alku == 0u)
							denominatorTOF(ax, local_ele, d_OSEM, local_ind, TOFSum, store_elements, DD[lor], TOFCenter, sigma_x, &D, ii * NBINS, d_epps, d_N);
#else
						if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && local_sino != 0.f && d_alku == 0u)
							denominator_multi(local_ele, &axCOSEM, &d_OSEM[local_ind]);
#endif
#else
						if (RHS) {
#ifdef TOF
							//denominatorTOF(ax, local_ele, d_OSEM, local_ind, TOFSum, store_elements, DD[lor], TOFCenter, sigma_x, &D, ii * NBINS, d_epps, d_N);
							denominatorTOF(ax, local_ele, d_OSEM, localInd, TOFSum, store_elements, DD[lor], TOFCenter, sigma_x, &D, ii* NBINS, d_epps, d_N);
#else
#ifdef AF
							//denominator(local_ele, ax, local_ind, d_N, d_OSEM);
							denominator(local_ele, ax, localInd, d_N, d_OSEM);
#else
							//denominator_multi(local_ele, &axOSEM, &d_OSEM[local_ind]);
							apuvar += read_imagef(d_OSEM, samplerIm, localInd).x;
							denominator_multi(local_ele, &axOSEM, read_imagef(d_OSEM, samplerIm, localInd).x);
#endif
#endif
						}
#endif
#endif
						//if (i.x == 100 && i.y == 100 && i.z == 0) {
						//	printf("axOSEM = %f\n", axOSEM);
						//		printf("tempi = %d\n", tempi);
						//		printf("tempj = %d\n", tempj);
						//		printf("tempk = %d\n", tempk);
						//}
						Np_n[lor]++;
						if (tempi < 0 || tempi >= d_Nxyz.x || tempj < 0 || tempj >= d_Nxyz.y || tempk < 0 || tempk >= d_Nxyz.z)
							break;
					}
					//pass[lor] = true;
					//if (i.x == 100 && i.y == 100 && i.z == 0) {
					//	printf("apuvar = %f\n", apuvar);
					//}
				}		
			}
			//if (i.x == 100 && i.y == 100 && i.z == 0) {
			//	printf("axOSEM = %f\n", axOSEM);
			//	printf("temp = %f\n", temp);
			//		//printf("tempi = %d\n", tempi);
			//		//printf("tempj = %d\n", tempj);
			//		//printf("tempk = %d\n", tempk);
			//}
		}
	}
	bool alku = true;
	for (lor = 0u; lor < N_RAYS; lor++) {
		if (pass[lor]) {
			if (alku) {
				//if (temp == 0.f)
				//	return;
				temp = 1.f / temp;
#ifdef ATN
				float n_r_summa = 0.f;
				for (ushort ln_r = 0u; ln_r < N_RAYS; ln_r++)
					n_r_summa += convert_float(pass[ln_r]);
				temp *= native_exp(jelppi / n_r_summa);
#endif
#ifdef NORM
				temp *= d_norm[idx];
#endif
#ifdef SCATTER
				temp *= d_scat[idx];
#endif
				temp *= global_factor;
#ifdef FP
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
					if (axCOSEM < d_epps)
						axCOSEM = d_epps;
					else
						axCOSEM *= temp;
#ifdef RANDOMS
					axCOSEM += d_sc_ra[idx];
#endif
					axCOSEM = local_sino / axCOSEM;
#endif
				}
#else
#ifdef AF
				if (RHS) {
#ifdef TOF
					nominatorTOF(MethodList, ax, d_Sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx, m_size, local_sino);
#else
					nominator(MethodList, ax, local_sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx);
#endif
				}
#else
				//if (i.x == 100 && i.y == 100 && i.z == 0) {
				//	printf("lor_toka = %d\n", lor);
				//	printf("axOSEM * temp = %f\n", axOSEM * temp);
				//}
				if (RHS) {
#ifdef TOF
#ifdef BP
					nominatorTOF(MethodList, ax, d_Sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx, m_size, local_sino);
#endif
					if (fp == 1) {
#pragma unroll NBINS
						for (int to = 0; to < NBINS; to++)
							d_output[idx + to * m_size + cumsum] = ax[to] * temp;
					}
#else
#ifdef BP
					if (axOSEM == 0.f) {
						axOSEM = d_epps;
					}
					else {
						axOSEM *= temp;
					}
#ifdef RANDOMS
					axOSEM += d_sc_ra[idx];
#endif
#endif
					if (fp == 1) {
						d_output[idx] = axOSEM * temp;
					}
					else {
						axOSEM = local_sino / axOSEM;
					}
#endif
				}
				if (fp == 1)
					return;
#endif
#endif
#endif
				alku = false;
			}
#ifdef BP
#ifdef TOF
			D = DD[lor];
#endif
			if (tx0_a[lor] > 1e6f && ty0_a[lor] > 1e6f) {
				const uint tempk = tempk_a[lor];
				//if (tempk >= d_N)
				//	return;
				if (ty0_a[lor] > tx0_a[lor]) {
					if (RHS) {
						for (uint k = 0; k < Np_n[lor]; k++) {
#ifdef TOF
#ifndef DEC
							const float TOFSum = TOFLoop(DD[lor], d_d.x, store_elements, TOFCenter, sigma_x, &D, k * NBINS, d_epps);
#endif
							backprojectTOF(tempk + k, d_d.x * temp, k * NBINS, store_elements, ax, d_Summ,
#ifndef DEC
								temp, sigma_x, & D, DD[lor], TOFCenter, d_epps, TOFSum,
#endif
#ifdef MBSREM
								MethodListOpenCL, d_alku, MBSREM_prepass, minimi, axACOSEM, d_OSEM, d_E, d_co, d_aco, local_sino, idx, m_size);
#else
								d_output, no_norm, d_N);
#endif
#else
#ifdef MBSREM

							if (d_alku == 0u) {
								if (MBSREM_prepass == 1)
#ifdef ATOMIC
									atom_add(&d_Summ[tempk + k], convert_long(d_d.x * temp * TH));
#elif defined(ATOMIC32)
									atomic_add(&d_Summ[tempk + k], convert_int(d_d.x * temp * TH));
#else
									atomicAdd_g_f(&d_Summ[tempk + k], (d_d.x * temp));
#endif
								if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
									minimi = d_d.x * temp;
									d_E[idx] += d_d.x * temp;
								}
								if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino != 0.f)
#ifdef ATOMIC
									atom_add(&d_co[tempk + k], convert_long(axCOSEM * d_d.x * temp * TH));
#elif defined(ATOMIC32)
									atomic_add(&d_co[tempk + k], convert_int(axCOSEM* d_d.x * temp * TH));
#else
									atomicAdd_g_f(&d_co[tempk + k], axCOSEM * d_d.x * temp);
#endif
								if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino != 0.f)
#ifdef ATOMIC
									atom_add(&d_aco[tempk + k], convert_long(axCOSEM * d_d.x * temp * TH));
#elif defined(ATOMIC32)
									atomic_add(&d_aco[tempk + k], convert_int(axCOSEM * d_d.x* temp * TH));
#else
									atomicAdd_g_f(&d_aco[tempk + k], axCOSEM * d_d.x * temp);
#endif
							}
							if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u)
								axACOSEM += (d_d.x * temp * read_imagef(d_OSEM, samplerIm, localInd).x);
#else
							if (no_norm == 0u)
#ifdef ATOMIC
								atom_add(&d_Summ[tempk + k], convert_long(d_d.x * temp * TH));
#elif defined(ATOMIC32)
								atomic_add(&d_Summ[tempk + k], convert_int(d_d.x* temp* TH));
#else
								atomicAdd_g_f(&d_Summ[tempk + k], (d_d.x * temp));
#endif
#if defined(FP) && defined(BP)
							if (local_sino != 0.f) {
#endif
#ifdef AF
								rhs(MethodList, d_d.x * temp, ax, tempk + k, d_N, d_output);
#else
#ifdef ATOMIC
								atom_add(&d_output[tempk + k], convert_long(d_d.x * temp * axOSEM * TH));
#elif defined(ATOMIC32)
								atomic_add(&d_output[tempk + k], convert_int(axOSEM * d_d.x * temp * TH));
#else
								atomicAdd_g_f(&d_output[tempk + k], ((d_d.x * temp) * axOSEM));
#endif
#endif
#if defined(FP) && defined(BP)
							}
#endif
#endif
#endif
						}
					}
//					else {
//						for (uint k = 0; k < Np_n[lor]; k++) {
//#ifdef TOF
//#ifndef DEC
//							const float TOFSum = TOFLoop(DD[lor], d_d.x, store_elements, TOFCenter, sigma_x, &D, k * NBINS, d_epps);
//#endif
//							sensTOF(tempk + k, d_d.x * temp, k * NBINS, store_elements, d_Summ,
//#ifndef DEC
//								temp, sigma_x, & D, DD[lor], TOFCenter, d_epps, TOFSum,
//#endif
//#ifdef MBSREM
//								MethodListOpenCL, d_alku, MBSREM_prepass, minimi, axACOSEM, d_OSEM, d_E, idx, m_size,
//#endif
//								no_norm);
//#else
//#ifdef MBSREM
//							if (d_alku == 0u && (MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
//								minimi = d_d.x * temp;
//								d_E[idx] += d_d.x * temp;
//							}
//#endif
//#ifdef ATOMIC
//							atom_add(&d_Summ[tempk + k], convert_long(d_d.x * temp * TH));
//#elif defined(ATOMIC32)
//							atomic_add(&d_Summ[tempk + k], convert_int(d_d.x * temp * TH));
//#else
//							atomicAdd_g_f(&d_Summ[tempk + k], (d_d.x * temp));
//#endif
//#endif
//						}
//					}
				}
				else {
					if (RHS) {
						for (uint k = 0; k < Np_n[lor]; k++) {
#ifdef TOF
#ifndef DEC
							const float TOFSum = TOFLoop(DD[lor], d_d.y, store_elements, TOFCenter, sigma_x, &D, k * NBINS, d_epps);
#endif
							backprojectTOF(tempk + k * d_N0, d_d.y * temp, k * NBINS, store_elements, ax, d_Summ,
#ifndef DEC
								temp, sigma_x, &D, DD[lor], TOFCenter, d_epps, TOFSum,
#endif
#ifdef MBSREM
								MethodListOpenCL, d_alku, MBSREM_prepass, minimi, axACOSEM, d_OSEM, d_E, d_co, d_aco, local_sino, idx, m_size);
#else
								d_output, no_norm, d_N);
#endif
#else
#ifdef MBSREM

							if (d_alku == 0u) {
								if (MBSREM_prepass == 1)
#ifdef ATOMIC
									atom_add(&d_Summ[tempk + k * d_N0], convert_long(d_d.y * temp * TH));
#elif defined(ATOMIC32)
									atomic_add(&d_Summ[tempk + k * d_N0], convert_int(d_d.y* temp* TH));
#else
									atomicAdd_g_f(&d_Summ[tempk + k * d_N0], (d_d.y * temp));
#endif
								if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
									minimi = d_d.y * temp;
									d_E[idx] += d_d.y * temp;
								}
								if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino != 0.f)
#ifdef ATOMIC
									atom_add(&d_co[tempk + k * d_N0], convert_long(axCOSEM * d_d.y * temp * TH));
#elif defined(ATOMIC32)
									atomic_add(&d_co[tempk + k * d_N0], convert_int(axCOSEM * d_d.y * temp * TH));
#else
									atomicAdd_g_f(&d_co[tempk + k * d_N0], axCOSEM * d_d.y * temp);
#endif
								if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino != 0.f)
#ifdef ATOMIC
									atom_add(&d_aco[tempk + k * d_N0], convert_long(axCOSEM * d_d.y * temp * TH));
#elif defined(ATOMIC32)
									atomic_add(&d_aco[tempk + k * d_N0], convert_int(axCOSEM* d_d.y* temp* TH));
#else
									atomicAdd_g_f(&d_aco[tempk + k * d_N0], axCOSEM * d_d.y * temp);
#endif
							}
							if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u)
								axACOSEM += (d_d.y * temp * read_imagef(d_OSEM, samplerIm, localInd).x);
#else
#if defined(FP) && defined(BP)
							if (local_sino != 0.f) {
#endif
#ifdef AF
								rhs(MethodList, d_d.y * temp, ax, tempk + k * d_N0, d_N, d_output);
#else
#ifdef ATOMIC
								atom_add(&d_output[tempk + k * d_N0], convert_long(d_d.y * temp * axOSEM * TH));
#elif defined(ATOMIC32)
								atomic_add(&d_output[tempk + k * d_N0], convert_int(axOSEM * d_d.y * temp * TH));
#else
								atomicAdd_g_f(&d_output[tempk + k * d_N0], ((d_d.y * temp) * axOSEM));
#endif
#endif
#if defined(FP) && defined(BP)
							}
#endif
							if (no_norm == 0u)
#ifdef ATOMIC
								atom_add(&d_Summ[tempk + k * d_N0], convert_long(d_d.y * temp * TH));
#elif defined(ATOMIC32)
								atomic_add(&d_Summ[tempk + k * d_N0], convert_int(d_d.y* temp* TH));
#else
								atomicAdd_g_f(&d_Summ[tempk + k * d_N0], (d_d.y * temp));
#endif
#endif
#endif
						}
					}
//					else {
//						for (uint k = 0; k < Np_n[lor]; k++) {
//#ifdef TOF
//#ifndef DEC
//							const float TOFSum = TOFLoop(DD[lor], d_d.y, store_elements, TOFCenter, sigma_x, &D, k * NBINS, d_epps);
//#endif
//							sensTOF(tempk + k * d_N0, d_d.y * temp, k * NBINS, store_elements, d_Summ,
//#ifndef DEC
//								temp, sigma_x, &D, DD[lor], TOFCenter, d_epps, TOFSum,
//#endif
//#ifdef MBSREM
//								MethodListOpenCL, d_alku, MBSREM_prepass, minimi, axACOSEM, d_OSEM, d_E, idx, m_size,
//#endif
//								no_norm);
//#else
//#ifdef MBSREM
//							if (d_alku == 0u && (MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
//								minimi = d_d.y * temp;
//								d_E[idx] += d_d.y * temp;
//							}
//#endif
//#ifdef ATOMIC
//							atom_add(&d_Summ[tempk + k * d_N0], convert_long(d_d.y * temp * TH));
//#elif defined(ATOMIC32)
//							atomic_add(&d_Summ[tempk + k * d_N0], convert_int(d_d.y* temp* TH));
//#else
//							atomicAdd_g_f(&d_Summ[tempk + k * d_N0], (d_d.y * temp));
//#endif
//#endif
//						}
//					}
				}
			}
			else {
				float tx0 = tx0_a[lor];
				float ty0 = ty0_a[lor];
				float tz0 = tz0_a[lor];
				const float txu = txu_a[lor];
				const float tyu = tyu_a[lor];
				const float tzu = tzu_a[lor];
				int tempi = tempi_a[lor];
				int tempj = tempj_a[lor];
				int tempk = tempk_a[lor];
				const int iu = iu_a[lor];
				const int ju = ju_a[lor];
				const int ku = ku_a[lor];
				float tc = tc_a[lor];
				float local_ele;
				if (RHS) {
					for (uint ii = 0u; ii < Np_n[lor]; ii++) {
						const uint local_ind = compute_ind(tempj, tempi, tempk, d_N0, d_N1, d_N, d_N0_a[lor], d_Nxy);
						if (tz0 < ty0 && tz0 < tx0) {
							local_ele = compute_element_2nd(&tz0, &tc, LL[lor], tzu, ku, &tempk, temp);
						}
						else if (ty0 < tx0 && ty0 <= tz0) {
							local_ele = compute_element_2nd(&ty0, &tc, LL[lor], tyu, ju, &tempj, temp);
						}
						else if (tx0 <= ty0 && tx0 <= tz0) {
							local_ele = compute_element_2nd(&tx0, &tc, LL[lor], txu, iu, &tempi, temp);
						}
#ifdef TOF
#ifndef DEC
						const float TOFSum = TOFLoop(DD[lor], local_ele / temp, store_elements, TOFCenter, sigma_x, &D, ii * NBINS, d_epps);
#endif
						backprojectTOF(local_ind, local_ele, ii* NBINS, store_elements, ax, d_Summ,
#ifndef DEC
							temp, sigma_x, & D, DD[lor], TOFCenter, d_epps, TOFSum,
#endif
#ifdef MBSREM
							MethodListOpenCL, d_alku, MBSREM_prepass, minimi, axACOSEM, d_OSEM, d_E, d_co, d_aco, local_sino, idx, m_size);
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
							axACOSEM += (local_ele * d_OSEM[local_ind]);
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
							atomic_add(&d_output[local_ind], convert_int(axOSEM * local_ele * TH));
#else
							atomicAdd_g_f(&d_output[local_ind], (local_ele * axOSEM));
#endif
#endif
#if defined(FP) && defined(BP)
						}
#endif
#endif
#endif
					}
				}
//				else {
//					for (uint ii = 0u; ii < Np_n[lor]; ii++) {
//						const uint local_ind = compute_ind(tempj, tempi, tempk, d_N0, d_N1, d_N, d_N0_a[lor], d_Nxy);
//						if (tz0 < ty0 && tz0 < tx0) {
//							local_ele = compute_element_2nd(&tz0, &tc, LL[lor], tzu, ku, &tempk, temp);
//						}
//						else if (ty0 < tx0 && ty0 <= tz0) {
//							local_ele = compute_element_2nd(&ty0, &tc, LL[lor], tyu, ju, &tempj, temp);
//						}
//						else if (tx0 <= ty0 && tx0 <= tz0) {
//							local_ele = compute_element_2nd(&tx0, &tc, LL[lor], txu, iu, &tempi, temp);
//						}
//#ifdef TOF
//#ifndef DEC
//						const float TOFSum = TOFLoop(DD[lor], local_ele / temp, store_elements, TOFCenter, sigma_x, &D, ii * NBINS, d_epps);
//#endif
//						sensTOF(local_ind, local_ele, ii* NBINS, store_elements, d_Summ,
//#ifndef DEC
//							temp, sigma_x, & D, DD[lor], TOFCenter, d_epps, TOFSum,
//#endif
//#ifdef MBSREM
//							MethodListOpenCL, d_alku, MBSREM_prepass, minimi, axACOSEM, d_OSEM, d_E, idx, m_size,
//#endif
//							no_norm);
//#else
//#ifdef MBSREM
//						if (d_alku == 0u) {
//							if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
//								if (local_ele < minimi && local_ele > 0.f)
//									minimi = local_ele;
//								d_E[idx] += local_ele;
//							}
//						}
//						if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u)
//							axACOSEM += (local_ele * d_OSEM[local_ind]);
//#endif
//#ifdef ATOMIC
//						atom_add(&d_Summ[local_ind], convert_long(local_ele * TH));
//#elif defined(ATOMIC32)
//						atomic_add(&d_Summ[local_ind], convert_int(local_ele* TH));
//#else
//						atomicAdd_g_f(&d_Summ[local_ind], local_ele);
//#endif
//#endif
//					}
//				}
			}
#endif
		}
	}
#ifdef MBSREM
	if (!alku) {
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
			d_ACOSEM_lhs[idx] = axACOSEM;
		}
#endif
	}
#endif
}


//#ifndef AF
//__kernel void summa(const __global CAST* d_Summ_device, __global CAST* d_Summ_local, const __global CAST* d_rhs_device, __global CAST* d_rhs_local,
//	const uint im_dim, const uchar no_norm) {
//
//	uint gid = get_global_id(0);
//
//	for (uint i = gid; i < im_dim; i += get_global_size(0)) {
//		if (no_norm == 0u)
//			d_Summ_local[i] += d_Summ_device[i];
//		d_rhs_local[i] += d_rhs_device[i];
//	}
//}
//
//__kernel void mlem(const __global CAST* d_Summ, const __global CAST* d_rhs, __global float* d_mlem, const uint im_dim, const float d_epps) {
//
//	uint gid = get_global_id(0);
//
//	for (uint i = gid; i < im_dim; i += get_global_size(0)) {
//#if defined(ATOMIC) || defined(ATOMIC32)
//		float rhs = convert_float(d_rhs[i]);
//		float Summ = convert_float(d_Summ[i]);
//		if (rhs != 0.f) {
//			if (Summ == 0.f)
//				d_mlem[i] = d_mlem[i] / d_epps * (rhs / TH);
//			else
//				d_mlem[i] = d_mlem[i] / (Summ / TH) * (rhs / TH);
//		}
//		else {
//			if (Summ != 0.f)
//				d_mlem[i] = d_mlem[i] / (Summ / TH) * d_epps;
//		}
//#else
//		if (d_rhs[i] != 0.f) {
//			if (d_Summ[i] == 0.f)
//				d_mlem[i] = d_mlem[i] / d_epps * d_rhs[i];
//			else
//				d_mlem[i] = d_mlem[i] / d_Summ[i] * d_rhs[i];
//		}
//		else {
//			if (d_Summ[i] != 0.f)
//				d_mlem[i] = d_mlem[i] / d_Summ[i] * d_epps;
//		}
//#endif
//	}
//}
//
//
//#ifdef PSF
//__kernel void Convolution3D(const __global CAST* input, __global CAST* output,
//	__constant float* convolution_window, int window_size_x, int window_size_y, int window_size_z) {
//	int4 ind = (int4)(get_global_id(0), get_global_id(1), get_global_id(2), 0);
//	int4 ind_uus = (0, 0, 0, 0);
//	float result = 0.f;
//	const uint Nyx = get_global_size(0) * get_global_size(1);
//	//int radius_x = floor((float)window_size_x / 2.0f);
//	//int radius_y = floor((float)window_size_y / 2.0f);
//	//int radius_z = floor((float)window_size_z / 2.0f);
//	int c = 0;
//	for (int k = -window_size_z; k <= window_size_z; k++) {
//		for (int j = -window_size_y; j <= window_size_y; j++) {
//			for (int i = -window_size_x; i <= window_size_x; i++) {
//				ind_uus.x = ind.x + i;
//				ind_uus.y = ind.y + j;
//				ind_uus.z = ind.z + k;
//				if (ind_uus.x >= get_global_size(0))
//					ind_uus.x = ind.x - i + 1;
//				else if (ind_uus.x < 0)
//					ind_uus.x = ind.x - (i + 1);
//				if (ind_uus.y >= get_global_size(1))
//					ind_uus.y = ind.y - j + 1;
//				else if (ind_uus.y < 0)
//					ind_uus.y = ind.y - (j + 1);
//				if (ind_uus.z >= get_global_size(2))
//					ind_uus.z = ind.z - k + 1;
//				else if (ind_uus.z < 0)
//					ind_uus.z = ind.z - (k + 1);
//				uint indeksi = ind_uus.x + ind_uus.y * get_global_size(0) + ind_uus.z * Nyx;
//#if defined(ATOMIC) || defined(ATOMIC32)
//				float p = convert_float(input[indeksi]) / TH;
//#else
//				float p = input[indeksi];
//#endif
//				p *= convolution_window[c];
//				result += p;
//				c++;
//			}
//		}
//	}
//#ifdef ATOMIC
//	output[ind.x + ind.y * get_global_size(0) + ind.z * Nyx] = convert_long(result * TH);
//#elif defined(ATOMIC32)
//	output[ind.x + ind.y * get_global_size(0) + ind.z * Nyx] = convert_int(result * TH);
//#else
//	output[ind.x + ind.y * get_global_size(0) + ind.z * Nyx] = result;
//#endif
//}
//
//
//__kernel void Convolution3D_f(const __global float* input, __global float* output,
//	__constant float* convolution_window, int window_size_x, int window_size_y, int window_size_z) {
//	int4 ind = (int4)(get_global_id(0), get_global_id(1), get_global_id(2), 0);
//	int4 ind_uus = (0, 0, 0, 0);
//	float result = 0.f;
//	const uint Nyx = get_global_size(0) * get_global_size(1);
//	//int radius_x = floor((float)window_size_x / 2.0f);
//	//int radius_y = floor((float)window_size_y / 2.0f);
//	//int radius_z = floor((float)window_size_z / 2.0f);
//	int c = 0;
//	for (int k = -window_size_z; k <= window_size_z; k++) {
//		for (int j = -window_size_y; j <= window_size_y; j++) {
//			for (int i = -window_size_x; i <= window_size_x; i++) {
//				ind_uus.x = ind.x + i;
//				ind_uus.y = ind.y + j;
//				ind_uus.z = ind.z + k;
//				if (ind_uus.x >= get_global_size(0))
//					ind_uus.x = ind.x - i + 1;
//				else if (ind_uus.x < 0)
//					ind_uus.x = ind.x - (i + 1);
//				if (ind_uus.y >= get_global_size(1))
//					ind_uus.y = ind.y - j + 1;
//				else if (ind_uus.y < 0)
//					ind_uus.y = ind.y - (j + 1);
//				if (ind_uus.z >= get_global_size(2))
//					ind_uus.z = ind.z - k + 1;
//				else if (ind_uus.z < 0)
//					ind_uus.z = ind.z - (k + 1);
//				uint indeksi = ind_uus.x + ind_uus.y * get_global_size(0) + ind_uus.z * Nyx;
//				float p = input[indeksi];
//				p *= convolution_window[c];
//				result += p;
//				c++;
//			}
//		}
//	}
//	output[ind.x + ind.y * get_global_size(0) + ind.z * Nyx] = result;
//}
//
//__kernel void vectorDiv(const __global float* input, __global float* output) {
//	uint id = get_global_id(0);
//	output[id] = output[id] / input[id];
//}
//
//__kernel void vectorMult(const __global float* input, __global float* output) {
//	uint id = get_global_id(0);
//	output[id] *= input[id];
//}
//#endif
//#endif
//
//#ifdef NLM_
//__kernel void NLM(__global float* grad, const __global float* u, const __global float* u_ref, __constant float* gaussian, const int search_window_x, const int search_window_y,
//	const int search_window_z, const int patch_window_x, const int patch_window_y, const int patch_window_z, const uint Nx, const uint Ny, const uint Nz,
//	const float h, const float epps, const int Nxy, const int min_x, const int max_x, const int min_y, const int max_y, const int min_z,
//	const int max_z, const int type) {
//
//	int n = get_global_id(0);
//	const int z = n / Nxy;
//	const int y = (n - z * Nxy) / convert_int(Nx);
//	const int x = n - z * Nxy - y * convert_int(Nx);
//	if (z < min_z || z >= max_z || x < min_x || x >= max_x || y < min_y || y >= max_y)
//		return;
//	float weight_sum = 0.f;
//	float output = 0.f;
//	const float uj = u[n];
//	for (int k = -search_window_z; k <= search_window_z; k++) {
//		const int z_n = z + k;
//		for (int j = -search_window_y; j <= search_window_y; j++) {
//			const int y_n = y + j;
//			for (int i = -search_window_x; i <= search_window_x; i++) {
//				const int x_n = x + i;
//				const int dim_n = z_n * Nxy + y_n * convert_int(Nx) + x_n;
//				const float uk = u[dim_n];
//				float distance = 0.f;
//				float weight = 0.f;
//
//				for (int pz = -patch_window_z; pz <= patch_window_z; pz++) {
//					const int z_k = (z_n + pz) * Nxy;
//					const int z_j = (z + pz) * Nxy;
//					for (int py = -patch_window_y; py <= patch_window_y; py++) {
//						const int y_k = (y_n + py) * convert_int(Nx);
//						const int y_j = (y + py) * convert_int(Nx);
//						int dim_g = (pz + patch_window_z) * (patch_window_x * 2 + 1) * (patch_window_y * 2 + 1) + (py + patch_window_y) * (patch_window_x * 2 + 1);
//						for (int px = -patch_window_x; px <= patch_window_x; px++) {
//							const float gg = gaussian[dim_g++];
//							//const float gg = 1.;
//							const int x_k = x_n + px;
//							const int dim_k = z_k + y_k + x_k;
//							const float Pj = u_ref[dim_k];
//							const int x_j = x + px;
//							const int dim_j = z_j + y_j + x_j;
//							const float Pk = u_ref[dim_j];
//							distance += gg * (Pj - Pk) * (Pj - Pk);
//						}
//					}
//				}
//				weight = exp(-distance / h);
//				weight_sum += weight;
//				if (type == 2)
//					output += weight * uk;
//				else if (type == 0) {
//					output += (weight * (uj - uk));
//				}
//				else {
//					output += ((weight * (uj - uk)) / sqrt(weight * (uj - uk) * (uj - uk) + epps));
//				}
//			}
//		}
//	}
//	weight_sum = 1.f / weight_sum;
//	output *= weight_sum;
//
//	grad[n] = output;
//
//}
//#endif
//
//#ifdef MEDIAN
//__kernel void medianFilter3D(const __global float* grad, __global float* output, const uint Nx, const uint Ny, const uint Nz) {
//	int xid = get_global_id(0);
//	int yid = get_global_id(1);
//	int zid = get_global_id(2);
//	if (xid < SEARCH_WINDOW_X || xid >= Nx + SEARCH_WINDOW_X || yid < SEARCH_WINDOW_Y || yid >= Ny + SEARCH_WINDOW_Y || zid < SEARCH_WINDOW_Z || zid >= Nz + SEARCH_WINDOW_Z)
//		return;
//	int koko = (SEARCH_WINDOW_X * 2 + 1) * (SEARCH_WINDOW_Y * 2 + 1) * (SEARCH_WINDOW_Z * 2 + 1);
//	float median[(SEARCH_WINDOW_X * 2 + 1) * (SEARCH_WINDOW_Y * 2 + 1) * (SEARCH_WINDOW_Z * 2 + 1)];
//	float medianF[(SEARCH_WINDOW_X * 2 + 1) * (SEARCH_WINDOW_Y * 2 + 1) * (SEARCH_WINDOW_Z * 2 + 1)];
//	int uu = 0;
//	for (int x = -SEARCH_WINDOW_X; x <= SEARCH_WINDOW_X; x++) {
//		for (int y = -SEARCH_WINDOW_Y; y <= SEARCH_WINDOW_Y; y++) {
//			for (int z = -SEARCH_WINDOW_Z; z <= SEARCH_WINDOW_Z; z++) {
//				uint pikseli = (xid + x) + (yid + y) * get_global_size(0) + (zid + z) * get_global_size(0) * get_global_size(1);
//				median[uu] = grad[pikseli];
//				uu++;
//			}
//		}
//	}
//	for (int hh = 0; hh < koko; hh++) {
//		int ind = 0;
//		for (int ll = 0; ll < koko; ll++) {
//			if (median[hh] > median[ll] || (median[hh] == median[ll] && hh < ll))
//				ind++;
//		}
//		medianF[ind] = median[hh];
//		if (ind == koko / 2)
//			break;
//	}
//	output[xid + yid * get_global_size(0) + zid * get_global_size(0) * get_global_size(1)] = medianF[koko / 2];
//}
//#endif