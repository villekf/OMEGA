/**************************************************************************
* Matrix free projectors for various ML and MAP algorithms.
* This function calculates the sensitivity image d_Summ = sum(A,1) (sum of 
* every row) and rhs = A*(y./(A'*x)) (forward and backward projections), 
* where A is the system matrix, y the measurements and x the 
* estimate/image.
* 
* This code is also used for the precomputation phase, where the number of
* voxels that each LOR traverses is computed.
*
* Used by implementations 2 and 3.
*
* This file contains all the three different projectors (Siddon, 
* orthogonal, volume-based). Preprocessor commands are used to separate
* different areas of code for the different projectors. This code also 
* includes the precomputation phase where the number of voxels in each LOR
* are computed. Furthermore the forward-backward projection example uses
* this same file. 64-bit atomics are also currently included in the same
* file and used if supported.
*
* Compiler preprocessing is utilized heavily, for example all the 
* corrections are implemented as compiler preprocesses. The code for 
* specific correction is thus only applied if it has been selected. The
* kernels are always compiled on-the-fly, though when using same input 
* parameters the kernel should be loaded from cache leading to a slightly
* faster startup time.
*
* INPUTS:
* global_factor = a global correction factor, e.g. dead time
* d_epps = a small constant to prevent division by zero,
* d_N = d_Nx * d_Ny * d_Nz,
* d_Nx/y/z = image size in x/y/z- dimension,
* d_dz/x/y = distance between adjecent voxels in z/x/y-dimension,
* d_bz/x/y = distance from the pixel space to origin (z/x/y-dimension),
* d_bzb = part in parenthesis of equation (9) in [1] precalculated when
* k = Nz,
* d_maxxx/yy = maximum distance of the pixel space from origin in
* x/y-dimension,
* d_zmax = maximum value of d_zdet,
* d_NSlices = the number of image slices,
* d_size_x = the number of detector elements,
* d_TotSinos = Total number of sinograms,
* d_det_per_ring = number of detectors per ring,
* d_raw = if 1 then raw data is used otherwise sinogram data
* pRows = number of pseudo rings,
* d_Nxy = d_Nx * d_Ny,
* fp = if 1, then only forward projection is computed, if 2 only 
* backprojection, if 0 then both,
* tube_width = the width of of the strip used for orthogonal distance based
* projector (2D),
* crystal_size_z = the width of of the tube used for orthogonal distance 
* based projector (3D),
* bmin = smaller orthogonal distances than this are fully inside the TOR, 
* volume projector only,
* bmax = Distances greater than this do not touch the TOR, volume projector
* only,
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
* x/y/z_center = Cartesian coordinates for the center of the voxels
* (x/y/z-axis),
* d_xy/zindex = for sinogram format they determine the detector
* indices corresponding to each sinogram bin (unused with raw data),
* V = Precomputed volumes for specific orthogonal distances, volume projector
* only
* d_L = detector numbers for raw data (unused for sinogram format),
* d_Sino = Sinogram/raw data,
* d_sc_ra = Randoms and/or scatter data,
* d_OSEM = buffer for OSEM/MLEM estimates,
* d_rhs_OSEM = buffer for OSEM/MLEM RHS elements,
* no_norm = If 1, normalization constant is not computed,
* m_size = Total number of LORs for this subset,
* cumsum = offset for input vector b in backprojection
*
* OUTPUTS:
* d_rhs_OSEM = RHS values for OSEM/MLEM,
* d_OSEM = OSEM/MLEM estimate,
* d_Summ = Sensitivity image
*
* [1] Jacobs, F., Sundermann, E., De Sutter, B., Christiaens, M. Lemahieu,
* I. (1998). A Fast Algorithm to Calculate the Exact Radiological Path
* through a Pixel or Voxel Space. Journal of computing and information
* technology, 6 (1), 89-94.
*
* Copyright (C) 2020 ViLe-Veikko Wettenhovi
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it wiL be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <https://www.gnu.org/licenses/>.
***************************************************************************/

// Matrix free orthogonal distance-based ray tracer, no precomputation step
__kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(LOCAL_SIZE, 1, 1)))
#ifdef FIND_LORS
void siddon_precomp(const uint d_Nxy, const uint d_N, const uint d_Nx, const uint d_Ny, const uint d_Nz, const float d_dz, const float d_dx, 
	const float d_dy, const float d_bz, const float d_bx, const float d_by, const float d_bzb, const float d_maxxx, const float d_maxyy, 
	const float d_zmax, const float d_NSlices, const uint d_size_x, const ushort d_TotSinos, const uint d_det_per_ring, const uchar d_raw, 
	const uint d_pRows, __constant uint* restrict d_pseudos, const __global float* restrict d_x, const __global float* restrict d_y, 
	const __global float* restrict d_zdet, __global ushort* restrict d_lor, const __global ushort* restrict d_L, const ulong m_size) {

#else

void kernel_multi(const float global_factor, const float d_epps, const uint d_N, const uint d_Nx, const uint d_Ny, const uint d_Nz, const float d_dz, 
	const float d_dx,	const float d_dy, const float d_bz, const float d_bx, const float d_by, const float d_bzb, const float d_maxxx, 
	const float d_maxyy, const float d_zmax, const float d_NSlices, const uint d_size_x, const uint d_TotSinos, const uint d_det_per_ring, 
	const uint d_pRows, const uint d_Nxy, const uchar fp, const float sigma_x, const float tube_width_xy, const float crystal_size_z, const float bmin, 
	const float bmax, const float Vmax, const float d_epsilon_mramla, __constant float* TOFCenter, const __global float* restrict d_atten, 
	__constant uint* d_pseudos, const __global float* restrict d_x, const __global float* restrict d_y, const __global float* restrict d_zdet,
	__constant float* x_center, __constant float* y_center, __constant float* z_center, __constant float* V, __constant uchar * MethodList,
#ifdef CT
	const uint d_subsets, __constant float* d_angles, const uint d_sizey, const float d_dPitch, const long d_nProjections, 
#endif
	const __global float* restrict d_norm, const __global float* restrict d_scat, __global CAST* restrict d_Summ, const __global ushort* restrict d_lor, 
	const __global uint* restrict d_xyindex, const __global ushort* restrict d_zindex, const __global ushort* restrict d_L, 
	const __global float* restrict d_Sino, const __global float* restrict d_sc_ra, const __global float* restrict d_OSEM,
#ifndef MBSREM
	__global CAST* restrict d_rhs_OSEM, const uchar no_norm, const ulong m_size, const ulong cumsum
#else
	const uint d_alku, const uchar MBSREM_prepass, __global float* restrict d_ACOSEM_lhs, __global float* restrict d_Amin, __global CAST* restrict d_co,
	__global CAST* restrict d_aco, __global float* restrict d_E, const ulong m_size, const RecMethodsOpenCL MethodListOpenCL, const ulong cumsum
#endif
) {
#endif
	// Get the current global index
	size_t idx = get_global_id(0);
	if (idx >= m_size)
		return;

#ifndef FIND_LORS // Not the precomputation phase
#ifdef TOF
	float local_sino = 0.f;
#ifndef LISTMODE2
#pragma unroll NBINS
	for (long to = 0L; to < NBINS; to++)
		local_sino += d_Sino[idx + m_size * to];
#endif
#else
#ifdef LISTMODE2
	const float local_sino = 0.f;
#else
	const float local_sino = (d_Sino[idx]);
#endif
#endif
#ifndef MBSREM
	if (no_norm == 1u && local_sino == 0.f)
		return;
#else
	const uchar no_norm = 0u;
#endif
#endif

	float xs, xd, ys, yd, zs, zd;
	// Load the next detector index
#ifdef CT
	get_detector_coordinates_CT(d_x, d_y, d_zdet, d_size_x, idx, d_subsets, d_angles, d_xyindex, d_zindex, d_sizey, d_dPitch, d_nProjections, 
		&xs, &xd, &ys, &yd, &zs, &zd);
#else
#ifdef RAW // raw data
#ifdef LISTMODE
	get_detector_coordinates(d_xyindex, d_zindex, d_size_x, idx, d_TotSinos, &xs, &xd, &ys, &yd, &zs, &zd, d_x, d_y, d_zdet, cumsum);
#else
	get_detector_coordinates_raw(d_x, d_y, d_zdet, d_L, d_det_per_ring, idx, d_pseudos, d_pRows, &xs, &xd, &ys, &yd, &zs, &zd);
#endif
#else // Sinogram data
#ifdef FIND_LORS // Precomputation phase
	get_detector_coordinates_precomp(d_size_x, idx, d_TotSinos, &xs, &xd, &ys, &yd, &zs, &zd, d_x, d_y, d_zdet);
#else // Not the precomputation phase
	get_detector_coordinates(d_xyindex, d_zindex, d_size_x, idx, d_TotSinos, &xs, &xd, &ys, &yd, &zs, &zd, d_x, d_y, d_zdet, cumsum);
#endif
#endif
#endif

	// Calculate the x, y and z distances of the detector pair
	float y_diff = (yd - ys);
	float x_diff = (xd - xs);
	const float z_diff = (zd - zs);

#ifdef PRECOMPUTE // Using precomputed data
	uint Np = convert_uint(d_lor[idx]);
#if !defined(DEC) || defined(TOF) // Intermediate results are not saved
	uint Np_n = Np;
#endif
#else // No precomputation
	if ((y_diff == 0.f && x_diff == 0.f && z_diff == 0.f) || (y_diff == 0.f && x_diff == 0.f))
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

#ifdef NORM // Normalization included
	local_norm = d_norm[idx];
#endif

	uint d_N0 = d_Nx;
	uint d_N1 = d_Ny;
	uint d_N2 = 1u;
	uint d_N3 = d_Nx;

#ifdef ORTH // 3D Orthogonal
	uint d_N4 = d_Nz;
#endif

	int tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 0;

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

#ifndef AF
#ifndef LISTMODE2
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

#ifdef ORTH // Orthogonal or volume-based ray tracer
	uchar xyz = 0u;
	float kerroin = 0.f;
	__constant float* xcenter = x_center;
	__constant float* ycenter = y_center;
#ifdef CRYST // 2.5D Orthogonal
	kerroin = e_norm(x_diff, y_diff, z_diff) * tube_width_xy;
#elif defined VOL // Volume-based
	kerroin = e_norm(x_diff, y_diff, z_diff);
#elif defined(CRYSTZ) && !defined(VOL) // 3D Orthogonal
	kerroin = e_norm(x_diff, y_diff, z_diff) * crystal_size_z;
#endif
#elif defined SIDDON // Siddon
	uint local_ind = 0u;
	float local_ele = 0.f;
#ifdef TOF
	float D = 0.f, DD = 0.f;
#endif
#endif
#else // Precomputation phase
	ushort temp_koko = 0u;
#endif
	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//If the LOR is perpendicular in the y-direction (Siddon cannot be used)
	if (fabs(z_diff) < 1e-6f && (fabs(y_diff) < 1e-6f || fabs(x_diff) < 1e-6f)) {

#ifdef FIND_LORS // Precomputation phase
		if (fabs(y_diff) < 1e-6f && yd <= d_maxyy && yd >= d_by && ys <= d_maxyy && ys >= d_by) {
			d_lor[idx] = convert_ushort(d_Nx);
			return;
		}
		else if (fabs(x_diff) < 1e-6f && xd <= d_maxxx && xd >= d_bx && xs <= d_maxxx && xs >= d_bx) {
			d_lor[idx] = convert_ushort(d_Ny);
			return;
		}
		else
			return;
#else // Not the precomputation phase

		tempk = convert_uint((zs - d_bz) / d_dz);

#ifdef SIDDON // Siddon

		float d_b, dd, d_d, d_d2;
		if (fabs(y_diff) < 1e-6f && yd <= d_maxyy && yd >= d_by && ys <= d_maxyy && ys >= d_by) {
			d_b = d_by;
			dd = yd;
			d_d = d_dy;
			d_d2 = d_dx;
			float xs_apu = xs;
			xs = ys;
			ys = xs_apu;
			float xdiff_apu = x_diff;
			x_diff = y_diff;
			y_diff = xdiff_apu;
			d_N0 = d_Ny;
			d_N1 = d_Nx;
			d_N2 = d_Ny;
			d_N3 = 1u;
		}
		else if (fabs(x_diff) < 1e-6f && xd <= d_maxxx && xd >= d_bx && xs <= d_maxxx && xs >= d_bx) {
			d_b = d_bx;
			dd = xd;
			d_d = d_dx;
			d_d2 = d_dy;
		}
		else
			return;
		float templ_ijk = 0.f;
		uint z_loop = 0u;
		perpendicular_elements(d_b, d_d, d_N0, dd, d_d2, d_N1, d_atten, &templ_ijk, &z_loop, tempk, d_N2, d_N3,
			d_norm, idx, global_factor, d_scat);
#ifdef TOF
		float dI = (d_d2 * d_N1) / 2.f * -sign(y_diff);
		D = dI;
		DD = D;
		float temp = templ_ijk / d_d2;
#ifdef DEC
		__private float store_elements[DEC * NBINS];
#else
		__private float store_elements[1];
#endif
#if defined(DEC) || defined(FP)
		local_ind = z_loop;
		for (uint ii = 0; ii < d_N1; ii++) {
			const float TOFSum = TOFLoop(DD, d_d2, store_elements, TOFCenter, sigma_x, &D, ii * NBINS, d_epps);
#ifdef FP 
			denominatorTOF(ax, d_d2, d_OSEM, local_ind, TOFSum, store_elements, DD, TOFCenter, sigma_x, &D, ii * NBINS, d_epps, d_N);
#endif
			local_ind += d_N3;
		}
#endif
#ifndef AF
#ifdef FP
		if (fp == 1) {
			nominatorTOF(MethodList, ax, d_Sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx, m_size, local_sino);
#pragma unroll NBINS
			for (long to = 0; to < NBINS; to++)
				d_rhs_OSEM[idx + to * m_size] = ax[to];
			return;
		}
#endif
#endif
#ifndef DEC
		D = DD;
#endif
		if (local_sino != 0.f) {
			local_ele = templ_ijk;
			local_ind = z_loop;
#ifdef FP
			nominatorTOF(MethodList, ax, d_Sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx, m_size, local_sino);
#endif
			for (uint ii = 0u; ii < d_N1; ii++) {
#ifndef DEC
				const float TOFSum = TOFLoop(DD, d_d2, store_elements, TOFCenter, sigma_x, &D, ii * NBINS, d_epps);
#endif
				backprojectTOF(local_ind, local_ele, ii * NBINS, store_elements, ax, d_Summ, 
#ifndef DEC
					temp, sigma_x, &D, DD, TOFCenter, d_epps, TOFSum,
#endif
#ifdef MBSREM
					MethodListOpenCL, d_alku, MBSREM_prepass, minimi, axACOSEM, d_OSEM, d_E, d_co, d_aco, local_sino, idx, m_size);
#else
					d_rhs_OSEM, no_norm, d_N);
#endif
				local_ind += d_N3;
			}
		}
		else {
			local_ele = templ_ijk;
			local_ind = z_loop;
			for (uint ii = 0u; ii < d_N1; ii++) {
#ifndef DEC
				const float TOFSum = TOFLoop(DD, d_d2, store_elements, TOFCenter, sigma_x, &D, ii * NBINS, d_epps);
#endif
				sensTOF(local_ind, local_ele, ii * NBINS, store_elements, d_Summ, 
#ifndef DEC
					temp, sigma_x, &D, DD, TOFCenter, d_epps, TOFSum,
#endif
#ifdef MBSREM
					MethodListOpenCL, d_alku, MBSREM_prepass, minimi, axACOSEM, d_OSEM, d_E, idx, m_size,
#endif
					no_norm);
				local_ind += d_N3;
			}
		}
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
#else
#ifdef MBSREM
		if (d_alku == 0u && ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0 ||
			(MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1 || MethodListOpenCL.RBIOSL == 1 || MethodListOpenCL.RBI == 1))) && local_sino != 0.f) {
			local_ele = templ_ijk;
			local_ind = z_loop;
			for (uint ii = 0u; ii < Np; ii++) {
				if (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0)
					axCOSEM += (local_ele * d_OSEM[local_ind]);
				if (MBSREM_prepass == 1)
#ifdef ATOMIC // 64-bit atomics
					atom_add(&d_Summ[local_ind], convert_long(local_ele * TH));
#elif defined(ATOMIC32)
					atomic_add(&d_Summ[local_ind], convert_int(local_ele * TH));
#else // 32-bit float atomics
					atomicAdd_g_f(&d_Summ[local_ind], local_ele);
#endif
				if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
					if (local_ele < minimi && local_ele > 0.f)
						minimi = local_ele;
					d_E[idx] += local_ele;
				}
				local_ind += d_N3;
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
			local_ind = z_loop;
			for (uint ii = 0u; ii < Np; ii++) {
				if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2))
#ifdef ATOMIC // 64-bit atomics
					atom_add(&d_co[local_ind], convert_long(axCOSEM * local_ele * TH));
#elif defined(ATOMIC32)
					atomic_add(&d_co[local_ind], convert_int(axCOSEM * local_ele * TH));
#else // 32-bit float atomics
					atomicAdd_g_f(&d_co[local_ind], axCOSEM * local_ele);
#endif
				if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1))
#ifdef ATOMIC // 64-bit atomics
					atom_add(&d_aco[local_ind], convert_long(axCOSEM * local_ele * TH));
#elif defined(ATOMIC32)
					atomic_add(&d_aco[local_ind], convert_int(axCOSEM * local_ele * TH));
#else // 32-bit float atomics
					atomicAdd_g_f(&d_aco[local_ind], axCOSEM * local_ele);
#endif
				local_ind += d_N3;
			}
		}
		else if (MBSREM_prepass == 1 && d_alku == 0u) {
			local_ele = templ_ijk;
			local_ind = z_loop;
			for (uint ii = 0u; ii < Np; ii++) {
#ifdef ATOMIC // 64-bit atomics
				atom_add(&d_Summ[local_ind], convert_long(local_ele * TH));
#elif defined(ATOMIC32)
				atomic_add(&d_Summ[local_ind], convert_int(local_ele* TH));
#else // 32-bit float atomics
				atomicAdd_g_f(&d_Summ[local_ind], local_ele);
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
			local_ind = z_loop;
			for (uint ii = 0u; ii < Np; ii++) {
				axACOSEM += (local_ele * d_OSEM[local_ind]);
				local_ind += d_N3;
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
#ifndef AF
		if (fp == 1) {
			local_ele = templ_ijk;
			local_ind = z_loop;
			for (uint ii = 0u; ii < d_N1; ii++) {
				denominator_multi(local_ele, &axOSEM, &d_OSEM[local_ind]);
				local_ind += d_N3;
			}
			nominator_multi(&axOSEM, local_sino, d_epps, 1.f, d_sc_ra, idx);
			d_rhs_OSEM[idx] = axOSEM;
			return;
		}
#endif
		// Calculate the next index and store it as weL as the probability of emission
		// If measurements are present, calculate the 
		if (local_sino != 0.f) {
			local_ele = templ_ijk;
			local_ind = z_loop;
#ifdef FP // Forward projection

			for (uint ii = 0u; ii < d_N1; ii++) {

#ifdef AF // Implementation 2

				denominator(local_ele, ax, local_ind, d_N, d_OSEM);

#else // Implementation 3

				denominator_multi(local_ele, &axOSEM, &d_OSEM[local_ind]);

#endif
				local_ind += d_N3;
			}
#ifdef AF // Implementation 2

			nominator(MethodList, ax, local_sino, d_epsilon_mramla, d_epps, 1.f, d_sc_ra, idx);

#else // Implementation 3

			nominator_multi(&axOSEM, local_sino, d_epps, 1.f, d_sc_ra, idx);

#endif
			local_ind = z_loop;

#endif
			
			for (uint ii = 0u; ii < d_N1; ii++) {
				if (local_ind >= 1032192)
					continue;
#ifdef AF
				rhs(MethodList, local_ele, ax, local_ind, d_N, d_rhs_OSEM);
#else
#ifdef ATOMIC // 64-bit atomics
				atom_add(&d_rhs_OSEM[local_ind], convert_long(local_ele * axOSEM * TH));
#elif defined(ATOMIC32)
				atomic_add(&d_rhs_OSEM[local_ind], convert_int(local_ele* axOSEM* TH));
#else // 32-bit float atomics
				atomicAdd_g_f(&d_rhs_OSEM[local_ind], (local_ele * axOSEM));
#endif
#endif
				if (no_norm == 0u)
#ifdef ATOMIC // 64-bit atomics
					atom_add(&d_Summ[local_ind], convert_long(local_ele * TH));
#elif defined(ATOMIC32)
					atomic_add(&d_Summ[local_ind], convert_int(local_ele* TH));
#else // 32-bit float atomics
					atomicAdd_g_f(&d_Summ[local_ind], local_ele);
#endif
				local_ind += d_N3;
			}
		}
		else {
			local_ele = templ_ijk;
			local_ind = z_loop;
			for (uint ii = 0u; ii < d_N1; ii++) {
#ifdef ATOMIC // 64-bit atomics
				atom_add(&d_Summ[local_ind], convert_long(local_ele * TH));
#elif defined(ATOMIC32)
				atomic_add(&d_Summ[local_ind], convert_int(local_ele* TH));
#else // 32-bit float atomics
				atomicAdd_g_f(&d_Summ[local_ind], local_ele);
#endif
				local_ind += d_N3;
			}
		}
#endif
#endif

#elif defined ORTH // Orthogonal or volume-based
		float temp = 0.f;
		float center2;
		float d_b, dd, d_d;
		if (fabs(y_diff) < 1e-6f && yd <= d_maxyy && yd >= d_by && ys <= d_maxyy && ys >= d_by) {
			center2 = x_center[0];
			d_b = d_by;
			dd = yd;
			d_d = d_dy;
			xcenter = y_center;
			float xs_apu = xs;
			xs = ys;
			ys = xs_apu;
			float xdiff_apu = x_diff;
			x_diff = y_diff;
			y_diff = xdiff_apu;
			d_N0 = d_Ny;
			d_N1 = d_Nx;
			d_N2 = d_Ny;
			d_N3 = 1u;
		}
		else if (fabs(x_diff) < 1e-6f && xd <= d_maxxx && xd >= d_bx && xs <= d_maxxx && xs >= d_bx) {
			center2 = y_center[0];
			d_b = d_bx;
			dd = xd;
			d_d = d_dx;
		}
		else
			return;

#ifdef CRYST // 2.5D orthogonal

#ifdef AF
#ifdef MBSREM
		orth_distance_perpendicular_multi(xcenter, center2, z_center, kerroin, &temp, &axACOSEM,
			d_b, dd, d_d, d_N0, d_N1, tempk, d_atten, local_norm, local_sino, xs, ys, zs, x_diff, y_diff, z_diff, d_N2, d_N3, d_OSEM, 0,
			d_Summ, true, false, global_factor, MethodListOpenCL, d_alku, &axCOSEM, d_E,
			d_co, d_aco, &minimi, MBSREM_prepass, d_sc_ra, d_Amin, d_ACOSEM_lhs, idx);
#else
		orth_distance_perpendicular_multi(xcenter, center2, z_center, kerroin, &temp, ax,
			d_b, dd, d_d, d_N0, d_N1, tempk, d_atten, local_norm, local_sino, xs, ys, zs, x_diff, y_diff, z_diff, d_N2, d_N3, d_OSEM, no_norm,
			d_Summ, true, false, global_factor, d_rhs_OSEM, d_N, MethodList);
#endif
#else
		orth_distance_perpendicular_multi(xcenter, center2, z_center, kerroin, &temp, &axOSEM,
			d_b, dd, d_d, d_N0, d_N1, tempk, d_atten, local_norm, local_sino, xs, ys, zs, x_diff, y_diff, z_diff, d_N2, d_N3, d_OSEM, no_norm,
			d_Summ, true, false, global_factor, d_rhs_OSEM, d_N, MethodList);

#endif

#else // 3D or volume-based

#ifdef AF
#ifdef MBSREM
		orth_distance_perpendicular_multi_3D(xcenter, center2, z_center, &temp, &axACOSEM,
			d_b, dd, d_d, d_N0, d_N1, tempk, d_atten, local_norm, local_sino, d_N2, d_N3, d_OSEM, xs, ys, zs, x_diff, y_diff, z_diff, kerroin, d_Nxy,
			d_Nz, 0, d_Summ, true, false, global_factor, bmin, bmax, Vmax, V, MethodListOpenCL, d_alku, &axCOSEM, d_E,
			d_co, d_aco, &minimi, MBSREM_prepass, d_sc_ra, d_Amin, d_ACOSEM_lhs, idx);
#else
		orth_distance_perpendicular_multi_3D(xcenter, center2, z_center, &temp, ax,
			d_b, dd, d_d, d_N0, d_N1, tempk, d_atten, local_norm, local_sino, d_N2, d_N3, d_OSEM, xs, ys, zs, x_diff, y_diff, z_diff, kerroin, d_Nxy,
			d_Nz, no_norm, d_Summ, true, false, global_factor, bmin, bmax, Vmax, V, d_rhs_OSEM, d_N, MethodList);
#endif
#else
		orth_distance_perpendicular_multi_3D(xcenter, center2, z_center, &temp, &axOSEM,
			d_b, dd, d_d, d_N0, d_N1, tempk, d_atten, local_norm, local_sino, d_N2, d_N3, d_OSEM, xs, ys, zs, x_diff, y_diff, z_diff, kerroin, d_Nxy,
			d_Nz, no_norm, d_Summ, true, false, global_factor, bmin, bmax, Vmax, V, d_rhs_OSEM, d_N, MethodList);
#endif

#endif

#ifndef AF
		if (fp == 1) {
			nominator_multi(&axOSEM, local_sino, d_epps, temp, d_sc_ra, idx);
			d_rhs_OSEM[idx] = axOSEM;
			return;
		}
#endif
#ifdef MBSREM
		if (d_alku == 0 && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1
			|| MethodListOpenCL.OSLCOSEM > 0) && local_sino != 0.f) {
			nominator_cosem(&axCOSEM, local_sino, d_epps, temp, d_sc_ra, idx);
		}
#ifdef CRYST // 2.5D orthogonal
		orth_distance_perpendicular_multi(xcenter, center2, z_center, kerroin, &temp, &axACOSEM,
			d_b, dd, d_d, d_N0, d_N1, tempk, d_atten, local_norm, local_sino, xs, ys, zs, x_diff, y_diff, z_diff, d_N2, d_N3, d_OSEM, 0,
			d_Summ, false, true, global_factor, MethodListOpenCL, d_alku, &axCOSEM, d_E,
			d_co, d_aco, &minimi, MBSREM_prepass, d_sc_ra, d_Amin, d_ACOSEM_lhs, idx);
#else
		orth_distance_perpendicular_multi_3D(xcenter, center2, z_center, &temp, &axACOSEM,
			d_b, dd, d_d, d_N0, d_N1, tempk, d_atten, local_norm, local_sino, d_N2, d_N3, d_OSEM, xs, ys, zs, x_diff, y_diff, z_diff, kerroin, d_Nxy,
			d_Nz, 0, d_Summ, false, true, global_factor, bmin, bmax, Vmax, V, MethodListOpenCL, d_alku, &axCOSEM, d_E,
			d_co, d_aco, &minimi, MBSREM_prepass, d_sc_ra, d_Amin, d_ACOSEM_lhs, idx);
#endif

#else
		if (local_sino != 0.f) {
#ifdef FP // Forward projection

#ifdef AF
			nominator(MethodList, ax, local_sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx);
#else
			nominator_multi(&axOSEM, local_sino, d_epps, temp, d_sc_ra, idx);
#endif

#endif

#ifdef CRYST // 2.5D orthogonal

#ifdef AF
			orth_distance_perpendicular_multi(xcenter, center2, z_center, kerroin, &temp, ax,
				d_b, dd, d_d, d_N0, d_N1, tempk, d_atten, local_norm, local_sino, xs, ys, zs, x_diff, y_diff, z_diff, d_N2, d_N3, d_OSEM, no_norm,
				d_Summ, false, true, global_factor, d_rhs_OSEM, d_N, MethodList);
#else
			orth_distance_perpendicular_multi(xcenter, center2, z_center, kerroin, &temp, &axOSEM,
				d_b, dd, d_d, d_N0, d_N1, tempk, d_atten, local_norm, local_sino, xs, ys, zs, x_diff, y_diff, z_diff, d_N2, d_N3, d_OSEM, no_norm,
				d_Summ, false, true, global_factor, d_rhs_OSEM, d_N, MethodList);
#endif

#else
#ifdef AF
			orth_distance_perpendicular_multi_3D(xcenter, center2, z_center, &temp, ax,
				d_b, dd, d_d, d_N0, d_N1, tempk, d_atten, local_norm, local_sino, d_N2, d_N3, d_OSEM, xs, ys, zs, x_diff, y_diff, z_diff, kerroin, d_Nxy,
				d_Nz, no_norm, d_Summ, false, true, global_factor, bmin, bmax, Vmax, V, d_rhs_OSEM, d_N, MethodList);
#else
			orth_distance_perpendicular_multi_3D(xcenter, center2, z_center, &temp, &axOSEM,
				d_b, dd, d_d, d_N0, d_N1, tempk, d_atten, local_norm, local_sino, d_N2, d_N3, d_OSEM, xs, ys, zs, x_diff, y_diff, z_diff, kerroin, d_Nxy,
				d_Nz, no_norm, d_Summ, false, true, global_factor, bmin, bmax, Vmax, V, d_rhs_OSEM, d_N, MethodList);
#endif
#endif
		}
		else {
#ifdef CRYST
#ifdef AF
			orth_distance_perpendicular_multi(xcenter, center2, z_center, kerroin, &temp, ax,
				d_b, dd, d_d, d_N0, d_N1, tempk, d_atten, local_norm, local_sino, xs, ys, zs, x_diff, y_diff, z_diff, d_N2, d_N3, d_OSEM, no_norm,
				d_Summ, false, false, global_factor, d_rhs_OSEM, d_N, MethodList);
#else
			orth_distance_perpendicular_multi(xcenter, center2, z_center, kerroin, &temp, &axOSEM,
				d_b, dd, d_d, d_N0, d_N1, tempk, d_atten, local_norm, local_sino, xs, ys, zs, x_diff, y_diff, z_diff, d_N2, d_N3, d_OSEM, no_norm,
				d_Summ, false, false, global_factor, d_rhs_OSEM, d_N, MethodList);
#endif

#else
#ifdef AF
			orth_distance_perpendicular_multi_3D(xcenter, center2, z_center, &temp, ax,
				d_b, dd, d_d, d_N0, d_N1, tempk, d_atten, local_norm, local_sino, d_N2, d_N3, d_OSEM, xs, ys, zs, x_diff, y_diff, z_diff, kerroin, d_Nxy,
				d_Nz, no_norm, d_Summ, false, false, global_factor, bmin, bmax, Vmax, V, d_rhs_OSEM, d_N, MethodList);
#else
			orth_distance_perpendicular_multi_3D(xcenter, center2, z_center, &temp, &axOSEM,
				d_b, dd, d_d, d_N0, d_N1, tempk, d_atten, local_norm, local_sino, d_N2, d_N3, d_OSEM, xs, ys, zs, x_diff, y_diff, z_diff, kerroin, d_Nxy,
				d_Nz, no_norm, d_Summ, false, false, global_factor, bmin, bmax, Vmax, V, d_rhs_OSEM, d_N, MethodList);
#endif
#endif
		}

#endif

#endif
#endif

	}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///*
	else {

		float txu = 0.f, tyu = 0.f, tzu = 0.f, tc = 0.f, tx0 = 1e8f, ty0 = 1e8f, tz0 = 1e8f;
		bool skip = false;

		// If the measurement is on a same ring
			// Z-coordinate (ring)
		if (fabs(z_diff) < 1e-6f) {
			tempk = convert_uint((zs - d_bz) / d_dz);
			skip = siddon_pre_loop_2D(d_bx, d_by, x_diff, y_diff, d_maxxx, d_maxyy, d_dx, d_dy, d_Nx, d_Ny, &tempi, &tempj, &txu, &tyu, &Np, TYPE,
				ys, xs, yd, xd, &tc, &iu, &ju, &tx0, &ty0);
		}
		else if (fabs(y_diff) < 1e-6f) {
			skip = siddon_pre_loop_2D(d_bx, d_bz, x_diff, z_diff, d_maxxx, d_bzb, d_dx, d_dz, d_Nx, d_Nz, &tempi, &tempk, &txu, &tzu, &Np, TYPE,
				zs, xs, zd, xd, &tc, &iu, &ku, &tx0, &tz0);
			tempj = perpendicular_start(d_by, yd, d_dy, d_Ny);
#ifndef PRECOMPUTE
			if (yd > d_maxyy || yd < d_by)
				skip = true;
#endif
		}
		else if (fabs(x_diff) < 1e-6f) {
			tempi = perpendicular_start(d_bx, xd, d_dx, d_Nx);
			skip = siddon_pre_loop_2D(d_by, d_bz, y_diff, z_diff, d_maxyy, d_bzb, d_dy, d_dz, d_Ny, d_Nz, &tempj, &tempk, &tyu, &tzu, &Np, TYPE,
				zs, ys, zd, yd, &tc, &ju, &ku, &ty0, &tz0);
#ifndef SIDDON
			int apu_tempi = tempi;
			float apu_txu = txu;
			float apu_tx0 = tx0;
			float apu_xdiff = x_diff;
			float apu_xs = xs;
			int apu_iu = iu;
			iu = ju;
			ju = apu_iu;
			tempi = tempj;
			tempj = apu_tempi;
			txu = tyu;
			tyu = apu_txu;
			tx0 = ty0;
			ty0 = apu_tx0;
			x_diff = y_diff;
			y_diff = apu_xdiff;
			xs = ys;
			ys = apu_xs;
			d_N0 = d_Ny;
			d_N1 = d_Nx;
			d_N2 = d_Ny;
			d_N3 = 1u;
#ifdef ORTH // Orthogonal or volume-based
			ycenter = x_center;
			xcenter = y_center;
#endif
#endif
#ifndef PRECOMPUTE
			if (xd > d_maxxx || xd < d_bx)
				skip = true;
#endif
		}
		else {
			skip = siddon_pre_loop_3D(d_bx, d_by, d_bz, x_diff, y_diff, z_diff, d_maxxx, d_maxyy, d_bzb, d_dx, d_dy, d_dz, d_Nx, d_Ny, d_Nz, &tempi, &tempj, &tempk, &tyu, &txu, &tzu,
				&Np, TYPE, ys, xs, yd, xd, zs, zd, &tc, &iu, &ju, &ku, &tx0, &ty0, &tz0);
		}
#ifndef PRECOMPUTE // No precomputation step performed
		if (skip)
			return;
#endif
#ifdef FIND_LORS // Precomputation phase
		for (uint ii = 0u; ii < Np; ii++) {
			temp_koko++;
			if (tz0 < ty0 && tz0 < tx0) {
				tempk += ku;
				tz0 += tzu;
			}
			else if (ty0 < tx0) {
				tempj += ju;
				ty0 += tyu;
			}
			else {
				tempi += iu;
				tx0 += txu;
			}
			if (tempj < 0 || tempi < 0 || tempk < 0 || tempi >= d_N0 || tempj >= d_N1 || tempk >= d_Nz)
				break;
		}
		d_lor[idx] = temp_koko;
		return;
#else // Not the precomputation phase
		float temp = 0.f;
#if defined(SIDDON) || !defined(DEC) // Siddon or no save of intermediate results
		const float tx0_a = tx0, ty0_a = ty0, tz0_a = tz0;
		const int tempi_a = tempi, tempj_a = tempj, tempk_a = tempk;
#endif
#if defined(SIDDON)// Siddon
		L = native_sqrt(x_diff * x_diff + z_diff * z_diff + y_diff * y_diff);
		const float tc_a = tc;
#ifdef TOF
		TOFDis(x_diff, y_diff, z_diff, tc, L, &D, &DD);
#endif
#if defined(DEC) && defined(TOF) // Save intermediate TOF results
		__private float store_elements[DEC * NBINS];
#elif defined(TOF)
		__private float store_elements[1];
#endif
#endif
#if defined(ATN) && !defined(SIDDON)
		L = native_sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);
#endif
#ifdef ORTH // Orthogonal or volume-based
		int tempi_b, tempj_b, tempk_b;
		uint ind = 0u;
#ifdef DEC // Save intermediate results
		__private float store_elements[DEC];
		__private uint store_indices[DEC];
#else
		__private float store_elements[1];
		__private uint store_indices[1];
#endif
#ifdef CRYST
		int alku = tempk + 1;
		int loppu = tempk;
#else
		int alku = d_Nz;
		int loppu = 0;
		if (ku > 0) {
			alku = tempk + 1;
		}
		else if (ku < 0) {
			loppu = tempk;
		}
#endif
		orth_distance_multi_3D(tempi, d_N0, d_N4, y_diff, x_diff, z_diff, y_center, x_center, z_center, &temp, d_N2, tempj, tempk,
			local_sino, d_OSEM, xs, ys, zs, d_Nxy, kerroin, no_norm, RHS, SUMMA, d_Summ, d_N1, d_N3, alku,
			iu, ju, loppu, bmin, bmax, Vmax, V, store_elements,
#ifdef AF
#ifdef MBSREM
			store_indices, &ind, &axACOSEM, MethodListOpenCL, d_E, idx, d_co, d_aco, &minimi, MBSREM_prepass, &axCOSEM, d_alku);
#else
			store_indices, & ind, d_rhs_OSEM, d_N, MethodList, ax);
#endif
#else
			store_indices, & ind, d_rhs_OSEM, d_N, MethodList, & axOSEM);
#endif
		for (uint ii = 0u; ii < Np; ii++) {
			if (tz0 < ty0 && tz0 < tx0) {
#ifdef ATN
				compute_attenuation(&tc, &jelppi, L, tz0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
#endif
				tempk += ku;
				if (tempk < d_Nz && tempk >= 0) {
					alku = tempk + 1;
					loppu = tempk;
					orth_distance_multi_3D(tempi, d_N0, d_N4, y_diff, x_diff, z_diff, y_center, x_center, z_center, &temp, d_N2, tempj, tempk,
						local_sino, d_OSEM, xs, ys, zs, d_Nxy, kerroin, no_norm, RHS, SUMMA, d_Summ, d_N1, d_N3, alku,
						iu, ju, loppu, bmin, bmax, Vmax, V, store_elements,
#ifdef AF
#ifdef MBSREM
						store_indices, &ind, &axACOSEM, MethodListOpenCL, d_E, idx, d_co, d_aco, &minimi, MBSREM_prepass, &axCOSEM, d_alku);
#else
						store_indices, & ind, d_rhs_OSEM, d_N, MethodList, ax);
#endif
#else
						store_indices, & ind, d_rhs_OSEM, d_N, MethodList, & axOSEM);
#endif
				}
				tz0 += tzu;
				xyz = 3u;
			}
			else if (ty0 < tx0) {
#ifdef ATN
				compute_attenuation(&tc, &jelppi, L, ty0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
#endif
				tempj += (ju);
				ty0 += tyu;
				xyz = 2u;
			}
			else {
#ifdef ATN
				compute_attenuation(&tc, &jelppi, L, tx0, tempi, tempj, tempk, d_Nx, d_Nxy, d_atten);
#endif
				tempi += iu;
				tx0 += txu;
				xyz = 1u;
			}
#ifndef PRECOMPUTE
#ifndef DEC
			Np_n++;
#endif
			if (tempj < 0 || tempi < 0 || tempk < 0 || tempi >= d_N0 || tempj >= d_N1 || tempk >= d_Nz) {
#ifdef CRYSTZ
				if (xyz < 3 && fabs(z_diff) >= 1e-6f) {
					if (xyz == 1)
						tempi -= iu;
					else if (xyz == 2)
						tempj -= ju;
					if ((tempk >= (d_Nz - 1) && ku > 0) || (tempk <= 0 && ku < 0))
						break;
					else
						tempk += ku;
					alku = d_Nz;
					loppu = 0;
					if (ku > 0) {
						loppu = tempk;
					}
					else if (ku < 0) {
						alku = tempk + 1;
					}
					orth_distance_multi_3D(tempi, d_N0, d_N4, y_diff, x_diff, z_diff, y_center, x_center, z_center, &temp, d_N2, tempj, tempk,
						local_sino, d_OSEM, xs, ys, zs, d_Nxy, kerroin, no_norm, RHS, SUMMA, d_Summ, d_N1, d_N3, alku,
						iu, ju, loppu, bmin, bmax, Vmax, V, store_elements,
#ifdef AF
#ifdef MBSREM
						store_indices, &ind, &axACOSEM, MethodListOpenCL, d_E, idx, d_co, d_aco, &minimi, MBSREM_prepass, &axCOSEM, d_alku);
#else
						store_indices, & ind, d_rhs_OSEM, d_N, MethodList, ax);
#endif
#else
						store_indices, & ind, d_rhs_OSEM, d_N, MethodList, & axOSEM);
#endif
				}
#endif
				break;
			}
#endif
		}
#ifdef PRECOMPUTE
#ifdef CRYSTZ
		if (xyz < 3 && fabs(z_diff) >= 1e-6f) {
			if (xyz == 1)
				tempi -= iu;
			else if (xyz == 2)
				tempj -= ju;
			if ((tempk >= (d_Nz - 1) && ku > 0) || (tempk <= 0 && ku < 0) || ku == 0) {}
			else {
				tempk += ku;
				alku = d_Nz;
				loppu = 0;
				if (ku > 0) {
					loppu = tempk;
				}
				else if (ku < 0) {
					alku = tempk + 1;
				}
				orth_distance_multi_3D(tempi, d_N0, d_N4, y_diff, x_diff, z_diff, y_center, x_center, z_center, &temp, d_N2, tempj, tempk,
					local_sino, d_OSEM, xs, ys, zs, d_Nxy, kerroin, no_norm, RHS, SUMMA, d_Summ, d_N1, d_N3, alku,
					iu, ju, loppu, bmin, bmax, Vmax, V, store_elements,
#ifdef AF
#ifdef MBSREM
					store_indices, &ind, &axACOSEM, MethodListOpenCL, d_E, idx, d_co, d_aco, &minimi, MBSREM_prepass, &axCOSEM, d_alku);
#else
					store_indices, & ind, d_rhs_OSEM, d_N, MethodList, ax);
#endif
#else
					store_indices, & ind, d_rhs_OSEM, d_N, MethodList, & axOSEM);
#endif
			}
		}
#endif
#endif
#elif defined(SIDDON)
#if !defined(CT) || (defined(FP) && defined(CT))
		for (uint ii = 0U; ii < Np; ii++) {
			local_ind = compute_ind(tempj, tempi, tempk, d_N0, d_N1, d_N, d_N3, d_Nxy);
			if (tz0 < ty0 && tz0 < tx0) {
				local_ele = compute_element(&tz0, &tc, L, tzu, ku, &tempk, &temp);
			}
			else if (ty0 < tx0) {
				local_ele = compute_element(&ty0, &tc, L, tyu, ju, &tempj, &temp);
			}
			else {
				local_ele = compute_element(&tx0, &tc, L, txu, iu, &tempi, &temp);
			}
#ifdef ATN
			jelppi += (local_ele * -d_atten[local_ind]);
#endif
#if defined(TOF) && (defined(DEC) || defined(FP))
			const float TOFSum = TOFLoop(DD, local_ele, store_elements, TOFCenter, sigma_x, &D, ii * NBINS, d_epps);
#endif
#ifdef FP
			if (local_sino != 0.f) {


#ifdef AF // Implementation 2

#ifdef MBSREM
#ifdef TOF
				if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0u)
					denominatorTOF(ax, local_ele, d_OSEM, local_ind, TOFSum, store_elements, DD, TOFCenter, sigma_x, &D, ii * NBINS, d_epps, d_N);
#else
				if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0u)
					denominator_multi(local_ele, &axCOSEM, &d_OSEM[local_ind]);
#endif
#else
#ifdef TOF
				denominatorTOF(ax, local_ele, d_OSEM, local_ind, TOFSum, store_elements, DD, TOFCenter, sigma_x, &D, ii * NBINS, d_epps, d_N);
#else
				denominator(local_ele, ax, local_ind, d_N, d_OSEM);
#endif

#endif

#else // Implementation 3
#ifdef TOF
				denominatorTOF(ax, local_ele, d_OSEM, local_ind, TOFSum, store_elements, DD, TOFCenter, sigma_x, &D, ii* NBINS, d_epps, d_N);
#else
				denominator_multi(local_ele, &axOSEM, &d_OSEM[local_ind]);
#endif

#endif
			}
#endif
#ifndef PRECOMPUTE
			Np_n++;
			if (tempj < 0 || tempi < 0 || tempk < 0 || tempi >= d_N0 || tempj >= d_N1 || tempk >= d_Nz) {
				break;
			}
#endif
		}
#else
		Np_n = Np;
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
#else
#if !defined(DEC)
		xyz = 0u;
#endif
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
#ifndef AF
#ifdef FP 
		if (fp == 1) {
#ifdef TOF
			nominatorTOF(MethodList, ax, d_Sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx, m_size, local_sino);
#pragma unroll NBINS
			for (int to = 0; to < NBINS; to++)
				d_rhs_OSEM[idx + to * m_size] = ax[to];
#else
			nominator_multi(&axOSEM, local_sino, d_epps, temp, d_sc_ra, idx);
			d_rhs_OSEM[idx] = axOSEM;
#endif
			return;
		}
#endif
#endif

		if (local_sino != 0.f) {
#ifdef FP
#ifdef TOF
			nominatorTOF(MethodList, ax, d_Sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx, m_size, local_sino);
#else
#ifdef AF
			nominator(MethodList, ax, local_sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx);
#else
			nominator_multi(&axOSEM, local_sino, d_epps, temp, d_sc_ra, idx);
#endif
#endif
#endif
			RHS = true;
		}
		else
			SUMMA = true;
#endif
		/* Add additional computations before backprojection here */
#ifdef ORTH
#ifdef DEC
		for (uint ii = 0u; ii < ind; ii++) {
#ifdef CT
			const float local_ele = store_elements[ii];
#else
			const float local_ele = store_elements[ii] * temp;
#endif
			const uint local_ind = store_indices[ii];
			if (RHS) {
#ifdef AF
				rhs(MethodList, local_ele, ax, local_ind, d_N, d_rhs_OSEM);
#else
#ifdef ATOMIC
				atom_add(&d_rhs_OSEM[local_ind], convert_long(local_ele * axOSEM * TH));
#elif defined(ATOMIC32)
				atomic_add(&d_rhs_OSEM[local_ind], convert_int(local_ele * axOSEM * TH));
#else
				atomicAdd_g_f(&d_rhs_OSEM[local_ind], (local_ele * axOSEM));
#endif
#endif
			}
			if (no_norm == 0u)
#ifdef ATOMIC
				atom_add(&d_Summ[local_ind], convert_long(local_ele* TH));
#elif defined(ATOMIC32)
				atomic_add(&d_Summ[local_ind], convert_int(local_ele * TH));
#else
				atomicAdd_g_f(&d_Summ[local_ind], local_ele);
#endif
		}
#else
#ifdef CRYST
		alku = tempk + 1;
		loppu = tempk;
#else
		alku = d_Nz;
		loppu = 0;
		if (ku > 0) {
			alku = tempk + 1;
		}
		else if (ku < 0) {
			loppu = tempk;
		}
#endif
		orth_distance_multi_3D(tempi, d_N0, d_N4, y_diff, x_diff, z_diff, y_center, x_center, z_center, &temp, d_N2, tempj, tempk,
			local_sino, d_OSEM, xs, ys, zs, d_Nxy, kerroin, no_norm, RHS, SUMMA, d_Summ, d_N1, d_N3, alku,
			iu, ju, loppu, bmin, bmax, Vmax, V, store_elements,
#ifdef AF
#ifdef MBSREM
			store_indices, &ind, &axACOSEM, MethodListOpenCL, d_E, idx, d_co, d_aco, &minimi, MBSREM_prepass, &axCOSEM, d_alku);
#else
			store_indices, & ind, d_rhs_OSEM, d_N, MethodList, ax);
#endif
#else
			store_indices, & ind, d_rhs_OSEM, d_N, MethodList, & axOSEM);
#endif
		for (uint ii = 0u; ii < Np_n; ii++) {
			if (tz0 < ty0 && tz0 < tx0) {
				tempk += ku;
				if (tempk < d_Nz && tempk >= 0) {
					alku = tempk + 1;
					loppu = tempk;
					orth_distance_multi_3D(tempi, d_N0, d_N4, y_diff, x_diff, z_diff, y_center, x_center, z_center, &temp, d_N2, tempj, tempk,
						local_sino, d_OSEM, xs, ys, zs, d_Nxy, kerroin, no_norm, RHS, SUMMA, d_Summ, d_N1, d_N3, alku,
						iu, ju, loppu, bmin, bmax, Vmax, V, store_elements,
#ifdef AF
#ifdef MBSREM
						store_indices, &ind, &axACOSEM, MethodListOpenCL, d_E, idx, d_co, d_aco, &minimi, MBSREM_prepass, &axCOSEM, d_alku);
#else
						store_indices, & ind, d_rhs_OSEM, d_N, MethodList, ax);
#endif
#else
						store_indices, & ind, d_rhs_OSEM, d_N, MethodList, & axOSEM);
#endif
				}
				tz0 += tzu;
				xyz = 3u;
			}
			else if (ty0 < tx0) {
				tempj += (ju);
				ty0 += tyu;
				xyz = 2u;
			}
			else {
				tempi += iu;
				tx0 += txu;
				xyz = 1u;
			}
		}
#ifdef CRYSTZ
		if (xyz < 3 && fabs(z_diff) >= 1e-6f) {
			if (xyz == 1)
				tempi -= iu;
			else if (xyz == 2)
				tempj -= ju;
			if ((tempk >= (d_Nz - 1) && ku > 0) || (tempk <= 0 && ku < 0)) {}
			else {
				tempk += ku;
				alku = d_Nz;
				loppu = 0;
				if (ku > 0) {
					loppu = tempk;
				}
				else if (ku < 0) {
					alku = tempk + 1;
				}
				orth_distance_multi_3D(tempi, d_N0, d_N4, y_diff, x_diff, z_diff, y_center, x_center, z_center, &temp, d_N2, tempj, tempk,
					local_sino, d_OSEM, xs, ys, zs, d_Nxy, kerroin, no_norm, RHS, SUMMA, d_Summ, d_N1, d_N3, alku,
					iu, ju, loppu, bmin, bmax, Vmax, V, store_elements,
#ifdef AF
#ifdef MBSREM
					store_indices, &ind, &axACOSEM, MethodListOpenCL, d_E, idx, d_co, d_aco, &minimi, MBSREM_prepass, &axCOSEM, d_alku);
#else
					store_indices, & ind, d_rhs_OSEM, d_N, MethodList, ax);
#endif
#else
					store_indices, & ind, d_rhs_OSEM, d_N, MethodList, & axOSEM);
#endif
			}
		}
#endif
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
		if (RHS) {
			for (uint ii = 0u; ii < Np_n; ii++) {
				local_ind = compute_ind(tempj, tempi, tempk, d_N0, d_N1, d_N, d_N3, d_Nxy);
				if (tz0 < ty0 && tz0 < tx0) {
					local_ele = compute_element_2nd(&tz0, &tc, L, tzu, ku, &tempk, temp);
				}
				else if (ty0 < tx0) {
					local_ele = compute_element_2nd(&ty0, &tc, L, tyu, ju, &tempj, temp);
				}
				else {
					local_ele = compute_element_2nd(&tx0, &tc, L, txu, iu, &tempi, temp);
				}
#ifdef TOF
#ifndef DEC
				const float TOFSum = TOFLoop(DD, local_ele / temp, store_elements, TOFCenter, sigma_x, &D, ii * NBINS, d_epps);
#endif
				backprojectTOF(local_ind, local_ele, ii * NBINS, store_elements, ax, d_Summ, 
#ifndef DEC
					temp, sigma_x, &D, DD, TOFCenter, d_epps, TOFSum,
#endif
#ifdef MBSREM
					MethodListOpenCL, d_alku, MBSREM_prepass, minimi, axACOSEM, d_OSEM, d_E, d_co, d_aco, local_sino, idx, m_size);
#else
					d_rhs_OSEM, no_norm, d_N);
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

#ifdef AF
				rhs(MethodList, local_ele, ax, local_ind, d_N, d_rhs_OSEM);
#else

#ifdef ATOMIC
				atom_add(&d_rhs_OSEM[local_ind], convert_long(local_ele * axOSEM * TH));
#elif defined(ATOMIC32)
				atomic_add(&d_rhs_OSEM[local_ind], convert_int(local_ele * axOSEM * TH));
#else
				atomicAdd_g_f(&d_rhs_OSEM[local_ind], (local_ele * axOSEM));
#endif
#endif
#endif
#endif
#if defined(CT) && !defined(FP)
				if (tempj < 0 || tempi < 0 || tempk < 0 || tempi >= d_N0 || tempj >= d_N1 || tempk >= d_Nz) {
					break;
				}
#endif
			}
		}
		else {
			for (uint ii = 0u; ii < Np_n; ii++) {
				local_ind = compute_ind(tempj, tempi, tempk, d_N0, d_N1, d_N, d_N3, d_Nxy);
				if (tz0 < ty0 && tz0 < tx0) {
					local_ele = compute_element_2nd(&tz0, &tc, L, tzu, ku, &tempk, temp);
				}
				else if (ty0 < tx0) {
					local_ele = compute_element_2nd(&ty0, &tc, L, tyu, ju, &tempj, temp);
				}
				else {
					local_ele = compute_element_2nd(&tx0, &tc, L, txu, iu, &tempi, temp);
				}
#ifdef TOF
#ifndef DEC
				const float TOFSum = TOFLoop(DD, local_ele / temp, store_elements, TOFCenter, sigma_x, &D, ii * NBINS, d_epps);
#endif
				sensTOF(local_ind, local_ele, ii* NBINS, store_elements, d_Summ, 
#ifndef DEC
					temp, sigma_x, & D, DD, TOFCenter, d_epps, TOFSum,
#endif
#ifdef MBSREM
					MethodListOpenCL, d_alku, MBSREM_prepass, minimi, axACOSEM, d_OSEM, d_E, idx, m_size,
#endif
					no_norm);
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
				}
				if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u)
					axACOSEM += (local_ele * d_OSEM[local_ind]);
#else
#ifdef ATOMIC
				atom_add(&d_Summ[local_ind], convert_long(local_ele* TH));
#elif defined(ATOMIC32)
				atomic_add(&d_Summ[local_ind], convert_int(local_ele* TH));
#else
				atomicAdd_g_f(&d_Summ[local_ind], local_ele);
#endif
#endif
#endif
#if defined(CT) && !defined(FP)
				if (tempj < 0 || tempi < 0 || tempk < 0 || tempi >= d_N0 || tempj >= d_N1 || tempk >= d_Nz) {
					break;
				}
#endif
			}
		}
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
//*/
	}
}


#ifndef AF
__kernel void summa(const __global CAST* d_Summ_device, __global CAST* d_Summ_local, const __global CAST* d_rhs_device, __global CAST* d_rhs_local,
	const uint im_dim, const uchar no_norm) {

	uint gid = get_global_id(0);

	for (uint i = gid; i < im_dim; i += get_global_size(0)) {
		if (no_norm == 0u)
			d_Summ_local[i] += d_Summ_device[i];
		d_rhs_local[i] += d_rhs_device[i];
	}
}


__kernel void mlem(const __global CAST* d_Summ, const __global CAST* d_rhs, __global float* d_mlem, const uint im_dim, const float d_epps) {

	uint gid = get_global_id(0);

	for (uint i = gid; i < im_dim; i += get_global_size(0)) {
#ifdef LISTMODE2
#if defined(ATOMIC) || defined(ATOMIC32)
		float Summ = convert_float(d_Summ[i]);
		d_mlem[i] = Summ / TH;
#else
		d_mlem[i] = d_Summ[i];
#endif
#else
#if defined(ATOMIC) || defined(ATOMIC32)
		float rhs = convert_float(d_rhs[i]);
		float Summ = convert_float(d_Summ[i]);
		if (rhs != 0.f) {
			if (Summ == 0.f)
				d_mlem[i] = d_mlem[i] / d_epps * (rhs / TH);
			else
				d_mlem[i] = d_mlem[i] / (Summ / TH) * (rhs / TH);
		}
		else {
			if (Summ != 0.f)
				d_mlem[i] = d_mlem[i] / (Summ / TH) * d_epps;
		}
#else
		if (d_rhs[i] != 0.f) {
			if (d_Summ[i] == 0.f)
				d_mlem[i] = d_mlem[i] / d_epps * d_rhs[i];
			else
				d_mlem[i] = d_mlem[i] / d_Summ[i] * d_rhs[i];
		}
		else {
			if (d_Summ[i] != 0.f)
				d_mlem[i] = d_mlem[i] / d_Summ[i] * d_epps;
		}
#endif
		if (d_mlem[i] < d_epps)
			d_mlem[i] = d_epps;
#endif
	}
}

#ifdef PSF
__kernel void Convolution3D(const __global CAST* input, __global CAST* output,
	__constant float* convolution_window, int window_size_x, int window_size_y, int window_size_z) {
	int4 ind = (int4)(get_global_id(0), get_global_id(1), get_global_id(2), 0);
	int4 ind_uus = (int4)(0, 0, 0, 0);
	float result = 0.f;
	const uint Nyx = get_global_size(0) * get_global_size(1);
	//int radius_x = floor((float)window_size_x / 2.0f);
	//int radius_y = floor((float)window_size_y / 2.0f);
	//int radius_z = floor((float)window_size_z / 2.0f);
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
			if (ind_uus.z >= get_global_size(2))
				ind_uus.z = get_global_size(2) - 1 - (ind_uus.z - get_global_size(2));
			//if (ind_uus.z < 0)
			//	ind_uus.z = ind.z - (k + 1);
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
				if (ind_uus.y >= get_global_size(1))
					ind_uus.y = get_global_size(1) - 1 - (ind_uus.y - get_global_size(1));
				//if (ind_uus.y < 0)
				//	ind_uus.y = ind.y - (j + 1);
			}
			ind_uus.y *= get_global_size(0);
			for (int i = (-window_size_x); i <= window_size_x; i++) {
				//int indx = convert_int(ind.x);
				//indx += i;
				if (ind.x < window_size_x) {
					if (i < -ind.x)
						ind_uus.x = abs(i) - 1 - ind.x;
					else
						ind_uus.x = i + ind.x;
				}
				else {
					ind_uus.x = ind.x + i;
					if (ind_uus.x >= get_global_size(0))
						ind_uus.x = get_global_size(0) - 1 - (ind_uus.x - get_global_size(0));
					//if (ind_uus.x < 0)
					//	ind_uus.x = ind.x - (i + 1);
				}
				//if (indx >= get_global_size(0))
				//	indx = ind.x - i + 1;
				//else if (indx < 0)
				//	indx = abs(convert_int(ind.x)) - 1;
				//int indeksi = indx + ind_uus.y + ind_uus.z;
				int indeksi = ind_uus.x + ind_uus.y + ind_uus.z;
#if defined(ATOMIC) || defined(ATOMIC32)
				float p = convert_float(input[indeksi]) / TH;
#else
				float p = input[indeksi];
#endif
				p *= convolution_window[c];
				result += p;
				c++;
			}
		}
	}
#ifdef ATOMIC
	output[ind.x + ind.y * get_global_size(0) + ind.z * Nyx] = convert_long(result * TH);
#elif defined(ATOMIC32)
	output[ind.x + ind.y * get_global_size(0) + ind.z * Nyx] = convert_int(result * TH);
#else
	output[ind.x + ind.y * get_global_size(0) + ind.z * Nyx] = result;
#endif
}


__kernel void Convolution3D_f(const __global float* input, __global float* output,
	__constant float* convolution_window, int window_size_x, int window_size_y, int window_size_z) {
	int4 ind = (int4)(get_global_id(0), get_global_id(1), get_global_id(2), 0);
	int4 ind_uus = (int4)(0, 0, 0, 0);
	float result = 0.f;
	const uint Nyx = get_global_size(0) * get_global_size(1);
	//int radius_x = floor((float)window_size_x / 2.0f);
	//int radius_y = floor((float)window_size_y / 2.0f);
	//int radius_z = floor((float)window_size_z / 2.0f);
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
			if (ind_uus.z >= get_global_size(2))
				ind_uus.z = get_global_size(2) - 1 - (ind_uus.z - get_global_size(2));
			//if (ind_uus.z < 0)
			//	ind_uus.z = ind.z - (k + 1);
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
				if (ind_uus.y >= get_global_size(1))
					ind_uus.y = get_global_size(1) - 1 - (ind_uus.y - get_global_size(1));
				//if (ind_uus.y < 0)
				//	ind_uus.y = ind.y - (j + 1);
			}
			ind_uus.y *= get_global_size(0);
			for (int i = (-window_size_x); i <= window_size_x; i++) {
				//int indx = convert_int(ind.x);
				//indx += i;
				if (ind.x < window_size_x) {
					if (i < -ind.x)
						ind_uus.x = abs(i) - 1 - ind.x;
					else
						ind_uus.x = i + ind.x;
				}
				else {
					ind_uus.x = ind.x + i;
					if (ind_uus.x >= get_global_size(0))
						ind_uus.x = get_global_size(0) - 1 - (ind_uus.x - get_global_size(0));
					//if (ind_uus.x < 0)
					//	ind_uus.x = ind.x - (i + 1);
				}
				//if (indx >= get_global_size(0))
				//	indx = ind.x - i + 1;
				//else if (indx < 0)
				//	indx = abs(convert_int(ind.x)) - 1;
				//int indeksi = indx + ind_uus.y + ind_uus.z;
				int indeksi = ind_uus.x + ind_uus.y + ind_uus.z;
				float p = input[indeksi];
				p *= convolution_window[c];
				result += p;
				c++;
			}
		}
	}
	output[ind.x + ind.y * get_global_size(0) + ind.z * Nyx] = result;
}

__kernel void vectorDiv(const __global float* input, __global float* output, const float epps) {
	uint id = get_global_id(0);
	output[id] = output[id] / (input[id] + epps);
}

__kernel void vectorMult(const __global float* input, __global float* output) {
	uint id = get_global_id(0);
	output[id] *= input[id];
}
#endif
#endif

#ifdef NLM_
__kernel void NLM(__global float* grad, const __global float* u, const __global float* u_ref, __constant float* gaussian, const int search_window_x, const int search_window_y,
	const int search_window_z, const int patch_window_x, const int patch_window_y, const int patch_window_z, const uint Nx, const uint Ny, const uint Nz, 
	const float h, const float epps, const int Nxy, const int min_x, const int max_x, const int min_y, const int max_y, const int min_z, 
	const int max_z, const int type) {

	int n = get_global_id(0);
	const int z = n / Nxy;
	const int y = (n - z * Nxy) / convert_int(Nx);
	const int x = n - z * Nxy - y * convert_int(Nx);
	if (z < min_z || z >= max_z || x < min_x || x >= max_x || y < min_y || y >= max_y)
		return;
	float weight_sum = 0.f;
	float output = 0.f;
	const float uj = u[n];
	for (int k = -search_window_z; k <= search_window_z; k++) {
		const int z_n = z + k;
		for (int j = -search_window_y; j <= search_window_y; j++) {
			const int y_n = y + j;
			for (int i = -search_window_x; i <= search_window_x; i++) {
				const int x_n = x + i;
				const int dim_n = z_n * Nxy + y_n * convert_int(Nx) + x_n;
				const float uk = u[dim_n];
				float distance = 0.f;
				float weight = 0.f;

				for (int pz = -patch_window_z; pz <= patch_window_z; pz++) {
					const int z_k = (z_n + pz) * Nxy;
					const int z_j = (z + pz) * Nxy;
					for (int py = -patch_window_y; py <= patch_window_y; py++) {
						const int y_k = (y_n + py) * convert_int(Nx);
						const int y_j = (y + py) * convert_int(Nx);
						int dim_g = (pz + patch_window_z) * (patch_window_x * 2 + 1) * (patch_window_y * 2 + 1) + (py + patch_window_y) * (patch_window_x * 2 + 1);
						for (int px = -patch_window_x; px <= patch_window_x; px++) {
							const float gg = gaussian[dim_g++];
							//const float gg = 1.;
							const int x_k = x_n + px;
							const int dim_k = z_k + y_k + x_k;
							const float Pj = u_ref[dim_k];
							const int x_j = x + px;
							const int dim_j = z_j + y_j + x_j;
							const float Pk = u_ref[dim_j];
							distance += gg * (Pj - Pk) * (Pj - Pk);
						}
					}
				}
				weight = exp(-distance / h);
				weight_sum += weight;
				if (type == 2)
					output += weight * uk;
				else if (type == 0) {
					output += (weight * (uj - uk));
				}
				else {
					output += ((weight * (uj - uk)) / sqrt(weight * (uj - uk) * (uj - uk) + epps));
				}
			}
		}
	}
	weight_sum = 1.f / weight_sum;
	output *= weight_sum;

	grad[n] = output;

}
#endif

#ifdef MEDIAN
__kernel void medianFilter3D(const __global float* grad, __global float* output, const uint Nx, const uint Ny, const uint Nz) {
	int xid = get_global_id(0);
	int yid = get_global_id(1);
	int zid = get_global_id(2);
	if (xid < SEARCH_WINDOW_X || xid >= Nx + SEARCH_WINDOW_X || yid < SEARCH_WINDOW_Y || yid >= Ny + SEARCH_WINDOW_Y || zid < SEARCH_WINDOW_Z || zid >= Nz + SEARCH_WINDOW_Z)
		return;
	int koko = (SEARCH_WINDOW_X * 2 + 1) * (SEARCH_WINDOW_Y * 2 + 1) * (SEARCH_WINDOW_Z * 2 + 1);
	float median[(SEARCH_WINDOW_X * 2 + 1) * (SEARCH_WINDOW_Y * 2 + 1) * (SEARCH_WINDOW_Z * 2 + 1)];
	float medianF[(SEARCH_WINDOW_X * 2 + 1) * (SEARCH_WINDOW_Y * 2 + 1) * (SEARCH_WINDOW_Z * 2 + 1)];
	for (int ll = 0; ll < koko; ll++) {
		medianF[ll] = 0.f;
		median[ll] = 0.f;
	}
	int uu = 0;
	for (int x = -SEARCH_WINDOW_X; x <= SEARCH_WINDOW_X; x++) {
		for (int y = -SEARCH_WINDOW_Y; y <= SEARCH_WINDOW_Y; y++) {
			for (int z = -SEARCH_WINDOW_Z; z <= SEARCH_WINDOW_Z; z++) {
				int pikseli = (xid + x) + (yid + y) * get_global_size(0) + (zid + z) * get_global_size(0) * get_global_size(1);
				median[uu] = grad[pikseli];
				uu++;
			}
		}
	}
	for (int hh = 0; hh < koko; hh++) {
		int ind = 0;
		for (int ll = 0; ll < koko; ll++) {
			if (median[hh] > median[ll] || (median[hh] == median[ll] && hh < ll))
				ind++;
		}
		medianF[ind] = median[hh];
		if (ind == koko / 2)
			break;
	}
	output[xid + yid * get_global_size(0) + zid * get_global_size(0) * get_global_size(1)] = medianF[koko / 2];
}
#endif