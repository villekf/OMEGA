/**************************************************************************
* Matrix free projectors (CUDA).
* This function calculates d_Summ = sum(A,1) (sum of every row) and
* rhs = A*(y./(A'*x)), where A is the system matrix, y the measurements
* and x the estimate/image.
*
* Used by implementation 2.
*
* This file contains all the three different projectors (Siddon,
* orthogonal, volume-based). Preprocessor commands are used to separate
* different areas of code for the different projectors. This code also
* includes the precomputation phase where the number of voxels in each LOR
* are computed. 64-bit atomics are also currently included in the same
* file and used if supported and selected.
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
* pRows = number of pseudo rings,
* d_Nxy = d_Nx * d_Ny,
* fp = if 1, then only forward projection is computed, if 2 only
* backprojection, if 0 then both
* tube_width = the width of of the strip used for orthogonal distance based
* projector (2D),
* crystal_size_z = the width of of the tube used for orthogonal distance based
* projector (3D),
* bmin = smaller orthogonal distances than this are fully inside the TOR,
* volume projector only
* bmax = Distances greater than this do not touch the TOR, volume projector
* only
* Vmax = Full volume of the spherical "voxel", volume projector only
* d_atten = attenuation data (images),
* d_norm = normalization coefficients,
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
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <https://www.gnu.org/licenses/>.
***************************************************************************/
#ifdef ATOMIC
#define TH 100000000000.f
#endif
#define THR 0.01f
#ifndef N_REKOS
#define N_REKOS 1
#endif
#ifndef NBINS
#define NBINS 1
#endif
#define NROLLS (N_REKOS * NBINS)
#ifdef PRECOMPUTE
#define TYPE 1
#else
#define TYPE 0
#endif
#include "general_cuda_functions.cuh"
#ifdef VOL
#define CC 1e3f
#endif
#ifdef PSF_LIMIT
#include "cuda_functions_psf.cuh"
#elif defined CRYST
#include "cuda_functions_orth25D.cuh"
#include "cuda_functions_orth3D.cuh"
#elif defined CRYSTZ
#include "cuda_functions_orth3D.cuh"
#endif

// Matrix free orthogonal distance-based ray tracer, no precomputation step
extern "C" __global__
void kernel_multi(const float global_factor, const float d_epps, const unsigned int d_N, const unsigned int d_Nx, const unsigned int d_Ny, const unsigned int d_Nz,
	const float d_dz, const float d_dx, const float d_dy, const float d_bz, const float d_bx, const float d_by, const float d_bzb, const float d_maxxx,
	const float d_maxyy, const float d_zmax, const float d_NSlices, const unsigned int d_size_x, const unsigned short d_TotSinos, 
	const unsigned int d_det_per_ring, const unsigned int d_pRows, const unsigned int d_Nxy, const unsigned char fp, const float sigma_x, 
	const float tube_width_xy, const float crystal_size_z, const float bmin, const float bmax, const float Vmax, const float d_epsilon_mramla,
	const float* TOFCenter, const float* d_atten, const unsigned int * d_pseudos, const float* d_x, const float* d_y, const float* d_zdet,
	const float* x_center, const float* y_center, const float* z_center, const float* V, const unsigned char * MethodList,
#ifdef CT
	const unsigned int d_subsets, const float* d_angles, const unsigned int d_sizey, const float d_dPitch, const long long int d_nProjections,
#endif
	const float* d_norm, const float* d_scat, CAST * d_Summ, const unsigned short * d_lor, const unsigned int * d_xyindex, const unsigned short * d_zindex,
	const unsigned short * d_L, const float* d_Sino, const float* d_sc_ra, const float* d_OSEM,
#ifndef MBSREM
	 CAST * d_rhs_OSEM, const unsigned char no_norm, const unsigned long long int m_size, const unsigned long long int cumsum
#else
	const unsigned int d_alku, const unsigned char MBSREM_prepass, float* d_ACOSEM_lhs, float* d_Amin, CAST * d_co,
	 CAST * d_aco, float* d_E, const unsigned long long int m_size, const RecMethodsOpenCL MethodListOpenCL, const unsigned long long int cumsum
#endif
	//void kernel_multi(float * d_Summ) {
) {
	// Get the current global index
	size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx >= m_size)
		return;

#ifdef TOF
	float local_sino = 0.f;
#ifndef LISTMODE2
#pragma unroll NBINS
	for (long long int to = 0LL; to < NBINS; to++)
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
	const unsigned char no_norm = 0u;
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

	get_detector_coordinates(d_xyindex, d_zindex, d_size_x, idx, d_TotSinos, &xs, &xd, &ys, &yd, &zs, &zd, d_x, d_y, d_zdet, cumsum);

#endif

#endif

	// Calculate the x, y and z distances of the detector pair
	float y_diff = (yd - ys);
	float x_diff = (xd - xs);
	const float z_diff = (zd - zs);

#ifdef PRECOMPUTE // Using precomputed data
	unsigned int Np = (unsigned int)(d_lor[idx]);

#ifndef DEC // Intermediate results are not saved

	unsigned int Np_n = Np;

#endif

#else // No precomputation

	if ((y_diff == 0.f && x_diff == 0.f && z_diff == 0.f) || (y_diff == 0.f && x_diff == 0.f))
		return;
	unsigned int Np = 0u;

#if !defined(DEC) || defined(TOF) // Intermediate results are not saved

	unsigned int Np_n = 0u;

#endif

#endif

	bool RHS = false, SUMMA = false;
	float L;

#ifdef ATN // Attenuation included

	float jelppi = 0.f;

#endif

	float local_norm = 0.f;
#ifdef NORM // Normalization included

	local_norm = d_norm[idx];

#endif

	unsigned int d_N0 = d_Nx;
	unsigned int d_N1 = d_Ny;
	unsigned int d_N2 = 1u;
	unsigned int d_N3 = d_Nx;

#ifdef ORTH // 3D Orthogonal
	unsigned int d_N4 = d_Nz;
#endif

	int tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 0;

#ifdef AF
#ifdef MBSREM
#ifdef TOF
	float axACOSEM[NBINS];
	float ax[NBINS];
#pragma unroll NBINS
	for (unsigned long long int to = 0LL; to < NBINS; to++) {
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
	for (unsigned long long int to = 0LL; to < NBINS; to++)
		minimi[to] = 1e8f;
#else
	float minimi = 1e8f;
#endif
#else
	float ax[NROLLS];
#pragma unroll
	for (unsigned int kk = 0; kk < NROLLS; kk++)
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

	unsigned char xyz = 0u;
	float kerroin = 0.f;
	const float* xcenter = x_center;
	const float* ycenter = y_center;

#ifdef CRYST // 2.5D Orthogonal

	kerroin = e_norm(x_diff, y_diff, z_diff) * tube_width_xy;

#elif defined VOL // Volume-based

	kerroin = e_norm(x_diff, y_diff, z_diff);

#elif defined(CRYSTZ) && !defined(VOL) // 3D Orthogonal

	kerroin = e_norm(x_diff, y_diff, z_diff) * crystal_size_z;

#endif

#elif defined SIDDON // Siddon

	unsigned int local_ind = 0u;
	float local_ele = 0.f;
#ifdef TOF
	float D = 0.f, DD = 0.f;
#endif

#endif
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//If the LOR is perpendicular in the y-direction (Siddon cannot be used)
	if (fabs(z_diff) < 1e-6f && (fabs(y_diff) < 1e-6f || fabs(x_diff) < 1e-6f)) {

		tempk = __float2int_rz((zs / d_zmax) * (d_NSlices - 1.f));

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
		unsigned int z_loop = 0u;
		perpendicular_elements(d_b, d_d, d_N0, dd, d_d2, d_N1, d_atten, &templ_ijk, &z_loop, tempk, d_N2, d_N3,
			d_norm, idx, global_factor);
#ifdef TOF
		float dI = (d_d2 * d_N1) / copysignf(-2.f, y_diff);
		D = dI;
		DD = D;
		float temp = templ_ijk / d_d2;
#ifdef DEC
		float store_elements[DEC * NBINS];
#else
		float store_elements[1];
#endif
		local_ind = z_loop;
		for (unsigned int ii = 0; ii < d_N1; ii++) {
			const float TOFSum = TOFLoop(DD, d_d2, store_elements, TOFCenter, sigma_x, &D, ii * NBINS, d_epps);
			denominatorTOF(ax, d_d2, d_OSEM, local_ind, TOFSum, store_elements, DD, TOFCenter, sigma_x, &D, ii * NBINS, d_epps, d_N);
			local_ind += d_N3;
		}
#ifndef DEC
		D = DD;
#endif
		if (local_sino != 0.f) {
			local_ele = templ_ijk;
			local_ind = z_loop;
			nominatorTOF(MethodList, ax, d_Sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx, m_size, local_sino);
			for (unsigned int ii = 0u; ii < d_N1; ii++) {
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
			for (unsigned int ii = 0u; ii < d_N1; ii++) {
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
		if (d_alku == 0u && (MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
#pragma unroll NBINS
			for (long long int to = 0LL; to < NBINS; to++) {
				d_Amin[idx + to * m_size] = minimi[to];
			}
		}
		else if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u) {
#pragma unroll NBINS
			for (long long int to = 0LL; to < NBINS; to++) {
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
			(MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1 || MethodListOpenCL.RBIOSL == 1))) && local_sino > 0.f) {
			local_ele = templ_ijk;
			local_ind = z_loop;
			for (unsigned int ii = 0u; ii < Np; ii++) {
				if (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0)
					axCOSEM += (local_ele * d_OSEM[local_ind]);
				if (MBSREM_prepass == 1)
#ifdef ATOMIC // 64-bit atomics
					atomicAdd(&d_Summ[local_ind], __float2ull_rn(local_ele * TH));
#else // 32-bit float atomics
					atomicAdd(&d_Summ[local_ind], local_ele);
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
					axCOSEM = expf(-axCOSEM) / local_sino;
#else
					axCOSEM = local_sino / axCOSEM;
#endif
			}
			local_ele = templ_ijk;
			local_ind = z_loop;
			for (unsigned int ii = 0u; ii < Np; ii++) {
				if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2))
#ifdef ATOMIC // 64-bit atomics
					atomicAdd(&d_co[local_ind], __float2ull_rn(axCOSEM * local_ele * TH));
#else // 32-bit float atomics
					atomicAdd(&d_co[local_ind], axCOSEM * local_ele);
#endif
				if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1))
#ifdef ATOMIC // 64-bit atomics
					atomicAdd(&d_aco[local_ind], __float2ull_rn(axCOSEM * local_ele * TH));
#else // 32-bit float atomics
					atomicAdd(&d_aco[local_ind], axCOSEM * local_ele);
#endif
				local_ind += d_N3;
			}
		}
		else if (MBSREM_prepass == 1 && d_alku == 0u) {
			local_ele = templ_ijk;
			local_ind = z_loop;
			for (unsigned int ii = 0u; ii < Np; ii++) {
#ifdef ATOMIC // 64-bit atomics
				atomicAdd(&d_Summ[local_ind], __float2ull_rn(local_ele * TH));
#else // 32-bit float atomics
				atomicAdd(&d_Summ[local_ind], local_ele);
#endif
				if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1)) {
					if (local_ele < minimi && local_ele > 0.f)
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
			for (unsigned int ii = 0u; ii < Np; ii++) {
				axACOSEM += (local_ele * d_OSEM[local_ind]);
				local_ind += d_N3;
			}
#ifdef RANDOMS
				axACOSEM += d_sc_ra[idx];
#endif
#ifdef CT
				d_ACOSEM_lhs[idx] = expf(-axACOSEM);
#else
				d_ACOSEM_lhs[idx] = axACOSEM;
#endif
		}
#else
		// Calculate the next index and store it as weL as the probability of emission
		// If measurements are present, calculate the 
		if (local_sino > 0.f) {
			local_ele = templ_ijk;
			local_ind = z_loop;

			for (unsigned int ii = 0u; ii < Np; ii++) {

				denominator(local_ele, ax, local_ind, d_N, d_OSEM);
				local_ind += d_N3;
			}

			nominator(MethodList, ax, local_sino, d_epsilon_mramla, d_epps, 1.f, d_sc_ra, idx);

			local_ind = z_loop;

			for (unsigned int ii = 0u; ii < Np; ii++) {
				rhs(MethodList, local_ele, ax, local_ind, d_N, d_rhs_OSEM);

				if (no_norm == 0u)
#ifdef ATOMIC // 64-bit atomics
					atomicAdd(&d_Summ[local_ind], __float2ull_rn(local_ele * TH));
#else // 32-bit float atomics
					atomicAdd(&d_Summ[local_ind], local_ele);
#endif
				local_ind += d_N3;
			}
		}
		else {
			local_ele = templ_ijk;
			local_ind = z_loop;
			for (unsigned int ii = 0u; ii < Np; ii++) {
#ifdef ATOMIC // 64-bit atomics
				atomicAdd(&d_Summ[local_ind], __float2ull_rn(local_ele * TH));
#else // 32-bit float atomics
				atomicAdd(&d_Summ[local_ind], local_ele);
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

#else // 3D or volume-based

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

#endif

#ifdef MBSREM
		if (d_alku == 0 && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1
			|| MethodListOpenCL.OSLCOSEM > 0) && local_sino > 0.f) {
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
		if (local_sino > 0.f) {
			nominator(MethodList, ax, local_sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx);

#ifdef CRYST // 2.5D orthogonal

			orth_distance_perpendicular_multi(xcenter, center2, z_center, kerroin, &temp, ax,
				d_b, dd, d_d, d_N0, d_N1, tempk, d_atten, local_norm, local_sino, xs, ys, zs, x_diff, y_diff, z_diff, d_N2, d_N3, d_OSEM, no_norm,
				d_Summ, false, true, global_factor, d_rhs_OSEM, d_N, MethodList);

#else
			orth_distance_perpendicular_multi_3D(xcenter, center2, z_center, &temp, ax,
				d_b, dd, d_d, d_N0, d_N1, tempk, d_atten, local_norm, local_sino, d_N2, d_N3, d_OSEM, xs, ys, zs, x_diff, y_diff, z_diff, kerroin, d_Nxy,
				d_Nz, no_norm, d_Summ, false, true, global_factor, bmin, bmax, Vmax, V, d_rhs_OSEM, d_N, MethodList);
#endif
		}
		else {
#ifdef CRYST
			orth_distance_perpendicular_multi(xcenter, center2, z_center, kerroin, &temp, ax,
				d_b, dd, d_d, d_N0, d_N1, tempk, d_atten, local_norm, local_sino, xs, ys, zs, x_diff, y_diff, z_diff, d_N2, d_N3, d_OSEM, no_norm,
				d_Summ, false, false, global_factor, d_rhs_OSEM, d_N, MethodList);

#else
			orth_distance_perpendicular_multi_3D(xcenter, center2, z_center, &temp, ax,
				d_b, dd, d_d, d_N0, d_N1, tempk, d_atten, local_norm, local_sino, d_N2, d_N3, d_OSEM, xs, ys, zs, x_diff, y_diff, z_diff, kerroin, d_Nxy,
				d_Nz, no_norm, d_Summ, false, false, global_factor, bmin, bmax, Vmax, V, d_rhs_OSEM, d_N, MethodList);
#endif
		}

#endif

#endif

	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//*
	else {

		float txu = 0.f, tyu = 0.f, tzu = 0.f, tc = 0.f, tx0 = 1e8f, ty0 = 1e8f, tz0 = 1e8f;
		bool skip = false;

		// If the measurement is on a same ring
			// Z-coordinate (ring)
		if (fabs(z_diff) < 1e-6f) {
			tempk = __float2int_rz((zs / d_zmax) * (d_NSlices - 1.f));
			skip = siddon_pre_loop_2D(d_bx, d_by, x_diff, y_diff, d_maxxx, d_maxyy, d_dx, d_dy, d_Nx, d_Ny, &tempi, &tempj, &txu, &tyu, &Np, TYPE,
				ys, xs, yd, xd, &tc, &iu, &ju, &tx0, &ty0);
			//return;
		}
		else if (fabs(y_diff) < 1e-6f) {
			tempj = perpendicular_start(d_by, yd, d_dy, d_Ny);
			skip = siddon_pre_loop_2D(d_bx, d_bz, x_diff, z_diff, d_maxxx, d_bzb, d_dx, d_dz, d_Nx, d_Nz, &tempi, &tempk, &txu, &tzu, &Np, TYPE,
				zs, xs, zd, xd, &tc, &iu, &ku, &tx0, &tz0);
			if (yd > d_maxyy || yd < d_by)
				skip = true;
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
			if (xd > d_maxxx || xd < d_bx)
				skip = true;
		}
		else {
			skip = siddon_pre_loop_3D(d_bx, d_by, d_bz, x_diff, y_diff, z_diff, d_maxxx, d_maxyy, d_bzb, d_dx, d_dy, d_dz, d_Nx, d_Ny, d_Nz, &tempi, &tempj, &tempk, &tyu, &txu, &tzu,
				&Np, TYPE, ys, xs, yd, xd, zs, zd, &tc, &iu, &ju, &ku, &tx0, &ty0, &tz0);
		}
#ifndef PRECOMPUTE // No precomputation step performed
		if (skip)
			return;
#endif
		float temp = 0.f;
#if defined(SIDDON) || !defined(DEC) // Siddon or no save of intermediate results
		const float tx0_a = tx0, ty0_a = ty0, tz0_a = tz0;
		const int tempi_a = tempi, tempj_a = tempj, tempk_a = tempk;
#endif
#if defined(SIDDON) || defined(PSF_LIMIT)// Siddon
		L = sqrtf(x_diff * x_diff + z_diff * z_diff + y_diff * y_diff);
		const float tc_a = tc;
#ifdef TOF
		TOFDis(x_diff, y_diff, z_diff, tc, L, &D, &DD);
#endif
#if defined(DEC) && defined(TOF) // Save intermediate TOF results
		float store_elements[DEC * NBINS];
#elif defined(TOF)
		float store_elements[1];
#endif
#endif
#if defined(ATN) && !defined(SIDDON) && !defined(PSF_LIMIT)
		L = sqrtf(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);
#endif
#ifdef ORTH // Orthogonal or volume-based
		unsigned int ind = 0u;
		int tempi_b, tempj_b, tempk_b;
#ifdef DEC // Save intermediate results
		float store_elements[DEC];
		unsigned int store_indices[DEC];
#else
		float store_elements[1];
		unsigned int store_indices[1];
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
#ifdef MBSREM
			store_indices, &ind, &axACOSEM, MethodListOpenCL, d_E, idx, d_co, d_aco, &minimi, MBSREM_prepass, &axCOSEM, d_alku);
#else
			store_indices, & ind, d_rhs_OSEM, d_N, MethodList, ax);
#endif
		for (unsigned int ii = 0u; ii < Np; ii++) {
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
#ifdef MBSREM
						store_indices, &ind, &axACOSEM, MethodListOpenCL, d_E, idx, d_co, d_aco, &minimi, MBSREM_prepass, &axCOSEM, d_alku);
#else
						store_indices, & ind, d_rhs_OSEM, d_N, MethodList, ax);
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
#ifdef MBSREM
						store_indices, &ind, &axACOSEM, MethodListOpenCL, d_E, idx, d_co, d_aco, &minimi, MBSREM_prepass, &axCOSEM, d_alku);
#else
						store_indices, & ind, d_rhs_OSEM, d_N, MethodList, ax);
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
				tempi = tempi + iu;
			else if (xyz == 2)
				tempj = tempj + ju;
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
#ifdef MBSREM
					store_indices, &ind, &axACOSEM, MethodListOpenCL, d_E, idx, d_co, d_aco, &minimi, MBSREM_prepass, &axCOSEM, d_alku);
#else
					store_indices, & ind, d_rhs_OSEM, d_N, MethodList, ax);
#endif
			}
		}
#endif
#endif
#elif defined(SIDDON)
//#if !defined(CT) || (defined(FP) && defined(CT))
		for (unsigned int ii = 0u; ii < Np; ii++) {
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
#ifdef TOF
			const float TOFSum = TOFLoop(DD, local_ele, store_elements, TOFCenter, sigma_x, &D, ii * NBINS, d_epps);
#endif
			if (local_sino > 0.f) {
#ifdef MBSREM
#ifdef TOF
				if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0u)
					denominatorTOF(ax, local_ele, d_OSEM, local_ind, TOFSum, store_elements, DD, TOFCenter, sigma_x, &D, ii * NBINS, d_epps, d_N);
#else
				if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0u)
					axCOSEM += (local_ele * d_OSEM[local_ind]);
#endif
#else
#ifdef TOF
				denominatorTOF(ax, local_ele, d_OSEM, local_ind, TOFSum, store_elements, DD, TOFCenter, sigma_x, &D, ii* NBINS, d_epps, d_N);
#else
				denominator(local_ele, ax, local_ind, d_N, d_OSEM);
#endif
#endif
			}
#ifndef PRECOMPUTE
			Np_n++;
			if (tempj < 0 || tempi < 0 || tempk < 0 || tempi >= d_N0 || tempj >= d_N1 || tempk >= d_Nz) {
				break;
			}
#endif
		}
//#else
//		Np_n = Np;
//#endif
#endif

#ifndef CT
		temp = 1.f / temp;
#ifdef ATN
		temp *= expf(jelppi);
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
#if defined(PSF_LIMIT)
		tc = tc_a;
#endif
#endif
#endif
#ifdef MBSREM
		if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && local_sino > 0.f && d_alku == 0u) {
#ifdef TOF
#pragma unroll NBINS
			for (long long int to = 0LL; to < NBINS; to++) {
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
			axCOSEM = expf(-axCOSEM) / local_sino;
#else
			axCOSEM = local_sino / axCOSEM;
#endif
#endif
		}
		RHS = true;
#else

		if (local_sino > 0.f) {
#ifdef TOF
			nominatorTOF(MethodList, ax, d_Sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx, m_size, local_sino);
#else
			nominator(MethodList, ax, local_sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx);
#endif
			RHS = true;
		}
		else
			SUMMA = true;
#endif
#ifdef ORTH
#ifdef DEC
		for (unsigned int ii = 0u; ii < ind; ii++) {
#ifdef CT
			const float local_ele = store_elements[ii];
#else
			const float local_ele = store_elements[ii] * temp;
#endif
			const unsigned int local_ind = store_indices[ii];
			if (RHS) {
				rhs(MethodList, local_ele, ax, local_ind, d_N, d_rhs_OSEM);
			}
			if (no_norm == 0u)
#ifdef ATOMIC
				atomicAdd(&d_Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
				atomicAdd(&d_Summ[local_ind], local_ele);
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
#ifdef MBSREM
			store_indices, &ind, &axACOSEM, MethodListOpenCL, d_E, idx, d_co, d_aco, &minimi, MBSREM_prepass, &axCOSEM, d_alku);
#else
			store_indices, & ind, d_rhs_OSEM, d_N, MethodList, ax);
#endif
		for (unsigned int ii = 0u; ii < Np_n; ii++) {
			if (tz0 < ty0 && tz0 < tx0) {
				tempk += ku;
				if (tempk < d_Nz && tempk >= 0) {
					alku = tempk + 1;
					loppu = tempk;
					orth_distance_multi_3D(tempi, d_N0, d_N4, y_diff, x_diff, z_diff, y_center, x_center, z_center, &temp, d_N2, tempj, tempk,
						local_sino, d_OSEM, xs, ys, zs, d_Nxy, kerroin, no_norm, RHS, SUMMA, d_Summ, d_N1, d_N3, alku,
						iu, ju, loppu, bmin, bmax, Vmax, V, store_elements,
#ifdef MBSREM
						store_indices, &ind, &axACOSEM, MethodListOpenCL, d_E, idx, d_co, d_aco, &minimi, MBSREM_prepass, &axCOSEM, d_alku);
#else
						store_indices, & ind, d_rhs_OSEM, d_N, MethodList, ax);
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
				tempi = tempi + iu;
			else if (xyz == 2)
				tempj = tempj + ju;
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
#ifdef MBSREM
					store_indices, &ind, &axACOSEM, MethodListOpenCL, d_E, idx, d_co, d_aco, &minimi, MBSREM_prepass, &axCOSEM, d_alku);
#else
					store_indices, & ind, d_rhs_OSEM, d_N, MethodList, ax);
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
			for (unsigned int ii = 0u; ii < Np_n; ii++) {
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
				backprojectTOF(local_ind, local_ele, ii* NBINS, store_elements, ax, d_Summ,
#ifndef DEC
					temp, sigma_x, & D, DD, TOFCenter, d_epps, TOFSum,
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
						atomicAdd(&d_Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
						atomicAdd(&d_Summ[local_ind], local_ele);
#endif
					if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
						if (local_ele < minimi && local_ele > 0.f)
							minimi = local_ele;
						d_E[idx] += local_ele;
					}
					if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
#ifdef ATOMIC
						atomicAdd(&d_co[local_ind], __float2ull_rn(axCOSEM * local_ele * TH));
#else
						atomicAdd(&d_co[local_ind], axCOSEM * local_ele);
#endif
					if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
#ifdef ATOMIC
						atomicAdd(&d_aco[local_ind], __float2ull_rn(axCOSEM * local_ele * TH));
#else
						atomicAdd(&d_aco[local_ind], axCOSEM * local_ele);
#endif
				}
				if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u)
					axACOSEM += (local_ele * d_OSEM[local_ind]);
#else
				if (no_norm == 0u)
#ifdef ATOMIC
					atomicAdd(&d_Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
					atomicAdd(&d_Summ[local_ind], local_ele);
#endif

				rhs(MethodList, local_ele, ax, local_ind, d_N, d_rhs_OSEM);

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
			for (unsigned int ii = 0u; ii < Np_n; ii++) {
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
						atomicAdd(&d_Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
						atomicAdd(&d_Summ[local_ind], local_ele);
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
				atomicAdd(&d_Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
				atomicAdd(&d_Summ[local_ind], local_ele);
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
		for (long long int to = 0LL; to < NBINS; to++) {
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
			d_ACOSEM_lhs[idx] = expf(-axACOSEM);
#else
			d_ACOSEM_lhs[idx] = axACOSEM;
#endif
		}
#endif
#endif
#endif
	}
	//*/
}


#ifdef NLM_
extern "C" __global__
void NLM(float* grad, const float* u, const float* u_ref, const float* gaussian, const int search_window_x, const int search_window_y,
	const int search_window_z, const int patch_window_x, const int patch_window_y, const int patch_window_z, const unsigned int Nx, const unsigned int Ny, const unsigned int Nz,
	const float h, const float epps, const int Nxy, const int min_x, const int max_x, const int min_y, const int max_y, const int min_z,
	const int max_z, const int type) {

	int n = threadIdx.x + blockIdx.x * blockDim.x;
	const int z = n / Nxy;
	const int y = (n - z * Nxy) / (int)(Nx);
	const int x = n - z * Nxy - y * (int)(Nx);
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
				const int dim_n = z_n * Nxy + y_n * (int)(Nx) + x_n;
				const float uk = u[dim_n];
				float distance = 0.f;
				float weight = 0.f;

				for (int pz = -patch_window_z; pz <= patch_window_z; pz++) {
					const int z_k = (z_n + pz) * Nxy;
					const int z_j = (z + pz) * Nxy;
					for (int py = -patch_window_y; py <= patch_window_y; py++) {
						const int y_k = (y_n + py) * (int)(Nx);
						const int y_j = (y + py) * (int)(Nx);
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
				weight = expf(-distance / h);
				weight_sum += weight;
				if (type == 2)
					output += weight * uk;
				else if (type == 0) {
					output += (weight * (uj - uk));
				}
				else {
					output += ((weight * (uj - uk)) / sqrtf(weight * (uj - uk) * (uj - uk) + epps));
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
extern "C" __global__
void medianFilter3D(const float* grad, float* output, const unsigned int Nx, const unsigned int Ny, const unsigned int Nz) {
	int xid = threadIdx.x + blockIdx.x * blockDim.x;
	int yid = threadIdx.y + blockIdx.y * blockDim.y;
	int zid = threadIdx.z + blockIdx.z * blockDim.z;
	unsigned int get_global_size[] = { gridDim.x * blockDim.x, gridDim.y * blockDim.y };
	if (xid < SEARCH_WINDOW_X || xid >= Nx + SEARCH_WINDOW_X || yid < SEARCH_WINDOW_Y || yid >= Ny + SEARCH_WINDOW_Y || zid < SEARCH_WINDOW_Z || zid >= Nz + SEARCH_WINDOW_Z)
		return;
	int koko = (SEARCH_WINDOW_X * 2 + 1) * (SEARCH_WINDOW_Y * 2 + 1) * (SEARCH_WINDOW_Z * 2 + 1);
	float median[(SEARCH_WINDOW_X * 2 + 1) * (SEARCH_WINDOW_Y * 2 + 1) * (SEARCH_WINDOW_Z * 2 + 1)];
	float medianF[(SEARCH_WINDOW_X * 2 + 1) * (SEARCH_WINDOW_Y * 2 + 1) * (SEARCH_WINDOW_Z * 2 + 1)];
	int uu = 0;
	for (int x = -SEARCH_WINDOW_X; x <= SEARCH_WINDOW_X; x++) {
		for (int y = -SEARCH_WINDOW_Y; y <= SEARCH_WINDOW_Y; y++) {
			for (int z = -SEARCH_WINDOW_Z; z <= SEARCH_WINDOW_Z; z++) {
				unsigned int pikseli = (xid + x) + (yid + y) * get_global_size[0] + (zid + z) * get_global_size[0] * get_global_size[1];
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
	output[xid + yid * get_global_size[0] + zid * get_global_size[0] * get_global_size[1]] = medianF[koko / 2];
	//output[xid + yid * get_global_size[0] + zid * get_global_size[0] * get_global_size[1]] = grad[xid + yid * get_global_size[0] + zid * get_global_size[0] * get_global_size[1]];
}
#endif