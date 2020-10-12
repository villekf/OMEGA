/**************************************************************************
* A matrix free improved Siddon's for multiple rays (CUDA).
* This function calculates Summ = sum(A,1) (sum of every row) and
* rhs = A*(y./(A'*x)), where A is the system matrix, y the measurements
* and x the estimate/image.
*
* Used by implementation 2.
*
* This version goes through all the LORs and determines on-the-fly if they
* intersect with the voxel space. Uses (optionally) multiple rays.
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
* d_rhs_OSEM = buffer for OSEM/MLEM RHS elements,
* no_norm = If 1, normalization constant is not computed,
* m_size = Total number of LORs for this subset.
*
* OUTPUTS:
* d_rhs_OSEM = rhs values for OSEM/MLEM,
* d_OSEM = OSEM/MLEM estimate,
* d_Summ = Sensitivity image
*
* [1] Jacobs, F., Sundermann, E., De Sutter, B., Christiaens, M. Lemahieu,
* I. (1998). A Fast Algorithm to Calculate the Exact Radiological Path
* through a Pixel or Voxel Space. Journal of computing and information
* technology, 6 (1), 89-94.
*
* Copyright (C) 2020 Ville-Veikko Wettenhovi
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <https://www.gnu.org/licenses/>.
***************************************************************************/
#ifdef ATOMIC
#define TH 100000000000.f
#endif
#include "general_cuda_functions.cuh"
#define TYPE 0
#ifndef N_REKOS
#define N_REKOS 1
#endif
#ifndef NBINS
#define NBINS 1
#endif
#define NROLLS (N_REKOS * NBINS)

// Matrix free Improved Siddon's algorithm
extern "C" __global__
void siddon_multi(const float global_factor, const float d_epps, const unsigned int d_N, const unsigned int d_Nx, const unsigned int d_Ny, const unsigned int d_Nz, const float d_dz, const float d_dx,
	const float d_dy, const float d_bz, const float d_bx, const float d_by, const float d_bzb, const float d_maxxx, const float d_maxyy,
	const float d_zmax, const float d_NSlices, const unsigned int d_size_x, const unsigned short d_TotSinos, 
	const unsigned int d_det_per_ring, const unsigned int d_pRows, const unsigned int d_Nxy, const unsigned char fp, const float sigma_x, const float dc_z, const unsigned short n_rays, const float d_epsilon_mramla,
	const float* TOFCenter, const float* d_atten, const unsigned int* d_pseudos, const float* d_x, const float* d_y, const float* d_zdet,
	const unsigned char* MethodList, const float* d_norm, const float* d_scat, CAST* d_Summ, const unsigned short* d_lor,
	const unsigned int* d_xyindex, const unsigned short* d_zindex, const unsigned short* d_L, const float* d_Sino, const float* d_sc_ra, const float* d_OSEM,
#ifndef MBSREM
	CAST* d_rhs_OSEM, const unsigned char no_norm, const unsigned long long int m_size
#else
	const unsigned int d_alku, const unsigned char MBSREM_prepass, float* d_ACOSEM_lhs, float* d_Amin, CAST* d_co,
	CAST* d_aco, float* d_E, const unsigned long long int m_size, const RecMethodsOpenCL MethodListOpenCL
#endif
) {
	// Get the current global index
	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx >= m_size)
		return;
#ifdef TOF
	float local_sino = 0.f;
#pragma unroll NBINS
	for (long long int to = 0LL; to < NBINS; to++)
		local_sino += d_Sino[idx + m_size * to];
#else
	const float local_sino = (d_Sino[idx]);
#endif
#ifndef MBSREM
	if (no_norm == 1u && local_sino == 0.f)
		return;
#else
	const unsigned char no_norm = 0u;
#endif

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
	bool RHS = true;
#else
	bool RHS = local_sino > 0.f ? true : false;
	float ax[NROLLS];
#pragma unroll
	for (unsigned int kk = 0; kk < NROLLS; kk++)
		ax[kk] = 0.f;
#endif
#ifdef TOF
	float D = 0.f;
#endif
	unsigned int d_N0 = d_Nx;
	unsigned int d_N1 = d_Ny;
	unsigned int d_N2 = 1u;
	unsigned int d_N3 = d_Nx;
	float jelppi = 0.f;
	float temp = 0.f;
	int tempi_a[N_RAYS];
	int tempj_a[N_RAYS];
	int tempk_a[N_RAYS];
	unsigned int d_N0_a[N_RAYS];
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
	unsigned int Np_n[N_RAYS];
	bool pass[N_RAYS];
#ifdef TOF
	float DD[N_RAYS];
#if defined(DEC) // Save intermediate TOF results
	float store_elements[DEC * NBINS];
#else
	float store_elements[1];
#endif
#endif
	// Load the next detector index
	// raw list-mode data
#pragma unroll N_RAYS
	for (unsigned short lor = 0u; lor < N_RAYS; lor++) {
		float xs, xd, ys, yd, zs, zd;
#ifdef RAW
		get_detector_coordinates_raw_multiray(d_x, d_y, d_zdet, d_L, d_det_per_ring, idx, d_pseudos, d_pRows, &xs, &xd, &ys, &yd, &zs, &zd, lor + 1u, dc_z);
		// Sinogram data
#else
		get_detector_coordinates_multiray(d_xyindex, d_zindex, d_size_x, idx, d_TotSinos, &xs, &xd, &ys, &yd, &zs, &zd, d_x, d_y, d_zdet, lor + 1u, dc_z);
#endif
		// Calculate the x, y and z distances of the detector pair
		float y_diff = (yd - ys);
		float x_diff = (xd - xs);
		const float z_diff = (zd - zs);
		if ((y_diff == 0.f && x_diff == 0.f && z_diff == 0.f) || (y_diff == 0.f && x_diff == 0.f)) {
			//pass[lor] = false;
			return;
		}
		pass[lor] = false;
		unsigned int Np = 0u;
		Np_n[lor] = 0u;
		// If the measurement is on a same ring
		if (fabs(z_diff) < 1e-6f && (fabs(y_diff) < 1e-6f || fabs(x_diff) < 1e-6f)) {
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
				tx0_a[lor] = 1e7f, ty0_a[lor] = 1e9f;
				pass[lor] = true;
			}
			else if (fabs(x_diff) < 1e-6f && xd <= d_maxxx && xd >= d_bx && xs <= d_maxxx && xs >= d_bx) {
				d_b = d_bx;
				dd = xd;
				d_d = d_dx;
				d_d2 = d_dy;
				tx0_a[lor] = 1e9f, ty0_a[lor] = 1e7f;
				pass[lor] = true;
			}
			if (pass[lor]) {
				Np_n[lor] = d_N1;
				unsigned int tempk = (unsigned int)((zs / d_zmax) * (d_NSlices - 1.f));
				unsigned int apu = 0u;
				const float element = perpendicular_elements_multiray(d_b, d_d, d_N0, dd, d_d2, d_N1, d_atten, &apu, tempk, d_N2, d_N3, &jelppi);
				temp += element;
				tempk_a[lor] = apu;
#ifdef TOF
				float dI = (d_d2 * d_N1) / copysignf(-2.f, y_diff);
				D = dI;
				DD[lor] = D;
				uint local_ind = apu;
				for (uint ii = 0; ii < d_N1; ii++) {
					const float TOFSum = TOFLoop(DD[lor], d_d2, store_elements, TOFCenter, sigma_x, &D, ii * NBINS, d_epps);
					denominatorTOF(ax, d_d2, d_OSEM, local_ind, TOFSum, store_elements, DD[lor], TOFCenter, sigma_x, &D, ii * NBINS, d_epps, d_N);
					local_ind += d_N3;
				}
#else
#ifdef MBSREM
				if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && local_sino > 0.f && d_alku == 0u) {
					for (unsigned int k = 0u; k < d_N1; k++)
						axCOSEM += (d_d * d_OSEM[apu + k * d_N3]);
				}
#else
				if (RHS) {
					for (unsigned int k = 0u; k < d_N1; k++) {
						denominator(d_d, ax, apu + k * d_N3, d_N, d_OSEM);
					}
				}
#endif
#endif
			}
		}
		else {
			int tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 0;
			float txu = 0.f, tyu = 0.f, tzu = 0.f, tc = 0.f, tx0 = 1e8f, ty0 = 1e8f, tz0 = 1e8f;
			bool skip = false;
			if (fabs(z_diff) < 1e-6f) {
				tempk = (int)((zs / d_zmax) * (d_NSlices - 1.f));
				skip = siddon_pre_loop_2D(d_bx, d_by, x_diff, y_diff, d_maxxx, d_maxyy, d_dx, d_dy, d_Nx, d_Ny, &tempi, &tempj, &txu, &tyu, &Np, TYPE,
					ys, xs, yd, xd, &tc, &iu, &ju, &tx0, &ty0);
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
				if (xd > d_maxxx || xd < d_bx)
					skip = true;
			}
			else {
				skip = siddon_pre_loop_3D(d_bx, d_by, d_bz, x_diff, y_diff, z_diff, d_maxxx, d_maxyy, d_bzb, d_dx, d_dy, d_dz, d_Nx, d_Ny, d_Nz, &tempi, &tempj, &tempk, &tyu, &txu, &tzu,
					&Np, TYPE, ys, xs, yd, xd, zs, zd, &tc, &iu, &ju, &ku, &tx0, &ty0, &tz0);
			}
			if (!skip) {
				pass[lor] = true;
				//return;
			}
			if (pass[lor]) {
				LL[lor] = sqrtf(x_diff * x_diff + z_diff * z_diff + y_diff * y_diff);
				tx0_a[lor] = tx0, ty0_a[lor] = ty0, tz0_a[lor] = tz0, tc_a[lor] = tc;
				txu_a[lor] = txu, tyu_a[lor] = tyu, tzu_a[lor] = tzu;
				tempi_a[lor] = tempi, tempj_a[lor] = tempj, tempk_a[lor] = tempk;
				d_N0_a[lor] = d_N0;
				//unsigned int temp_ijk = compute_ind(tempj, tempi, tempk, d_N0, d_N1, d_N, d_N3, d_Nxy);
				//temp_ijk_a[lor] = temp_ijk;
				iu_a[lor] = iu, ju_a[lor] = ju, ku_a[lor] = ku;
#ifdef TOF
				TOFDis(x_diff, y_diff, z_diff, tc, LL[lor], &D, &DD[lor]);
#endif
				float local_ele;
				for (unsigned int ii = 0u; ii < Np; ii++) {
					const unsigned int local_ind = compute_ind(tempj, tempi, tempk, d_N0, d_N1, d_N, d_N0, d_Nxy);
					if (tz0 < ty0 && tz0 < tx0) {
						local_ele = compute_element(&tz0, &tc, LL[lor], tzu, ku, &tempk, &temp);
					}
					else if (ty0 < tx0) {
						local_ele = compute_element(&ty0, &tc, LL[lor], tyu, ju, &tempj, &temp);
					}
					else {
						local_ele = compute_element(&tx0, &tc, LL[lor], txu, iu, &tempi, &temp);
					}
#ifdef ATN
					jelppi += (local_ele * -d_atten[local_ind]);
#endif
#ifdef TOF
					const float TOFSum = TOFLoop(DD[lor], local_ele, store_elements, TOFCenter, sigma_x, &D, ii * NBINS, d_epps);
#endif
#ifdef MBSREM
#ifdef TOF
					if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && local_sino != 0.f && d_alku == 0u)
						denominatorTOF(ax, local_ele, d_OSEM, local_ind, TOFSum, store_elements, DD[lor], TOFCenter, sigma_x, &D, ii * NBINS, d_epps, d_N);
#else
					if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && local_sino != 0.f && d_alku == 0u)
						axCOSEM += (local_ele * d_OSEM[local_ind]);
#endif
#else
					if (RHS) {
#ifdef TOF
						denominatorTOF(ax, local_ele, d_OSEM, local_ind, TOFSum, store_elements, DD[lor], TOFCenter, sigma_x, &D, ii * NBINS, d_epps, d_N);
#else
						denominator(local_ele, ax, local_ind, d_N, d_OSEM);
#endif
					}
#endif
					Np_n[lor]++;
					if (tempj < 0 || tempi < 0 || tempk < 0 || tempi >= d_Nx || tempj >= d_Ny || tempk >= d_Nz)
						break;
				}
				//pass[lor] = true;
			}
		}
	}
	bool alku = true;
	for (unsigned short lor = 0u; lor < N_RAYS; lor++) {
		if (pass[lor]) {
			if (alku) {
				temp = 1.f / temp;
#ifdef ATN
				float n_r_summa = 0.f;
				for (unsigned short ln_r = 0u; ln_r < N_RAYS; ln_r++)
					n_r_summa += (float)(pass[ln_r]);
				temp *= expf(jelppi / n_r_summa);
#endif
#ifdef NORM
				temp *= d_norm[idx];
#endif
				temp *= global_factor;
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
				if (RHS) {
#ifdef TOF
					nominatorTOF(MethodList, ax, d_Sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx, m_size, local_sino);
#else
					nominator(MethodList, ax, local_sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx);
#endif
				}
#endif
				alku = false;
			}
#ifdef TOF
			D = DD[lor];
#endif
			if (tx0_a[lor] > 1e6f && ty0_a[lor] > 1e6f) {
				const unsigned int tempk = tempk_a[lor];
				if (ty0_a[lor] > tx0_a[lor]) {
					if (RHS) {
						for (unsigned int k = 0; k < Np_n[lor]; k++) {
#ifdef TOF
#ifndef DEC
							const float TOFSum = TOFLoop(DD[lor], d_dx, store_elements, TOFCenter, sigma_x, &D, k * NBINS, d_epps);
#endif
							backprojectTOF(tempk + k, d_dx * temp, k * NBINS, store_elements, ax, d_Summ,
#ifndef DEC
								temp, sigma_x, &D, DD[lor], TOFCenter, d_epps, TOFSum,
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
									atomicAdd(&d_Summ[tempk + k], __float2ull_rn(d_dx * temp * TH));
#else
									atomicAdd(&d_Summ[tempk + k], (d_dx * temp));
#endif
								if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
									minimi = d_dx * temp;
									d_E[idx] += d_dx * temp;
								}
								if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
#ifdef ATOMIC
									atomicAdd(&d_co[tempk + k], __float2ull_rn(axCOSEM * d_dx * temp * TH));
#else
									atomicAdd(&d_co[tempk + k], axCOSEM * d_dx * temp);
#endif
								if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
#ifdef ATOMIC
									atomicAdd(&d_aco[tempk + k], __float2ull_rn(axCOSEM * d_dx * temp * TH));
#else
									atomicAdd(&d_aco[tempk + k], axCOSEM * d_dx * temp);
#endif
							}
							if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u)
								axACOSEM += (d_dx * temp * d_OSEM[tempk + k]);
#else
							if (no_norm == 0u)
#ifdef ATOMIC
								atomicAdd(&d_Summ[tempk + k], __float2ull_rn(d_dx * temp * TH));
#else
								atomicAdd(&d_Summ[tempk + k], (d_dx * temp));
#endif
							rhs(MethodList, d_dx * temp, ax, tempk + k, d_N, d_rhs_OSEM);
#endif
#endif
						}
					}
					else {
						for (unsigned int k = 0; k < Np_n[lor]; k++) {
#ifdef TOF
#ifndef DEC
							const float TOFSum = TOFLoop(DD[lor], d_dx, store_elements, TOFCenter, sigma_x, &D, k * NBINS, d_epps);
#endif
							sensTOF(tempk + k, d_dx * temp, k * NBINS, store_elements, d_Summ,
#ifndef DEC
								temp, sigma_x, &D, DD[lor], TOFCenter, d_epps, TOFSum,
#endif
#ifdef MBSREM
								MethodListOpenCL, d_alku, MBSREM_prepass, minimi, axACOSEM, d_OSEM, d_E, idx, m_size,
#endif
								no_norm);
#else
#ifdef MBSREM
							if (d_alku == 0u && (MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
								minimi = d_dx * temp;
								d_E[idx] += d_dx * temp;
							}
#endif
#ifdef ATOMIC
							atomicAdd(&d_Summ[tempk + k], __float2ull_rn(d_dx * temp * TH));
#else
							atomicAdd(&d_Summ[tempk + k], (d_dx * temp));
#endif
#endif
						}
					}
				}
				else {
					if (RHS) {
						for (unsigned int k = 0; k < Np_n[lor]; k++) {
#ifdef TOF
#ifndef DEC
							const float TOFSum = TOFLoop(DD[lor], d_dy, store_elements, TOFCenter, sigma_x, &D, k * NBINS, d_epps);
#endif
							backprojectTOF(tempk + k * d_N0, d_dy * temp, k * NBINS, store_elements, ax, d_Summ,
#ifndef DEC
								temp, sigma_x, &D, DD[lor], TOFCenter, d_epps, TOFSum,
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
									atomicAdd(&d_Summ[tempk + k * d_N0], __float2ull_rn(d_dy * temp * TH));
#else
									atomicAdd(&d_Summ[tempk + k * d_N0], (d_dy * temp));
#endif
								if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
									minimi = d_dy * temp;
									d_E[idx] += d_dy * temp;
								}
								if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
#ifdef ATOMIC
									atomicAdd(&d_co[tempk + k * d_N0], __float2ull_rn(axCOSEM * d_dy * temp * TH));
#else
									atomicAdd(&d_co[tempk + k * d_N0], axCOSEM * d_dy * temp);
#endif
								if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
#ifdef ATOMIC
									atomicAdd(&d_aco[tempk + k * d_N0], __float2ull_rn(axCOSEM * d_dy * temp * TH));
#else
									atomicAdd(&d_aco[tempk + k * d_N0], axCOSEM * d_dy * temp);
#endif
							}
							if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u)
								axACOSEM += (d_dy * temp * d_OSEM[tempk + k * d_N0]);
#else
							rhs(MethodList, d_dy * temp, ax, tempk + k * d_N0, d_N, d_rhs_OSEM);
							if (no_norm == 0u)
#ifdef ATOMIC
								atomicAdd(&d_Summ[tempk + k * d_N0], __float2ull_rn(d_dy * temp * TH));
#else
								atomicAdd(&d_Summ[tempk + k * d_N0], (d_dy * temp));
#endif
#endif
#endif
						}
					}
					else {
						for (unsigned int k = 0; k < Np_n[lor]; k++) {
#ifdef TOF
#ifndef DEC
							const float TOFSum = TOFLoop(DD[lor], d_dy, store_elements, TOFCenter, sigma_x, &D, k * NBINS, d_epps);
#endif
							sensTOF(tempk + k * d_N0, d_dy * temp, k * NBINS, store_elements, d_Summ,
#ifndef DEC
								temp, sigma_x, &D, DD[lor], TOFCenter, d_epps, TOFSum,
#endif
#ifdef MBSREM
								MethodListOpenCL, d_alku, MBSREM_prepass, minimi, axACOSEM, d_OSEM, d_E, idx, m_size,
#endif
								no_norm);
#else
#ifdef MBSREM
							if (d_alku == 0u && (MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
								minimi = d_dy * temp;
								d_E[idx] += d_dy * temp;
							}
#endif
#ifdef ATOMIC
							atomicAdd(&d_Summ[tempk + k * d_N0], __float2ull_rn(d_dy * temp * TH));
#else
							atomicAdd(&d_Summ[tempk + k * d_N0], (d_dy * temp));
#endif
#endif
						}
					}
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
					for (unsigned int ii = 0u; ii < Np_n[lor]; ii++) {
						const unsigned int local_ind = compute_ind(tempj, tempi, tempk, d_N0, d_N1, d_N, d_N0_a[lor], d_Nxy);
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
					}
				}
				else {
					for (unsigned int ii = 0u; ii < Np_n[lor]; ii++) {
						const unsigned int local_ind = compute_ind(tempj, tempi, tempk, d_N0, d_N1, d_N, d_N0_a[lor], d_Nxy);
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
						sensTOF(local_ind, local_ele, ii* NBINS, store_elements, d_Summ,
#ifndef DEC
							temp, sigma_x, & D, DD[lor], TOFCenter, d_epps, TOFSum,
#endif
#ifdef MBSREM
							MethodListOpenCL, d_alku, MBSREM_prepass, minimi, axACOSEM, d_OSEM, d_E, idx, m_size,
#endif
							no_norm);
#else
#ifdef MBSREM
						if (d_alku == 0u) {
							if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
								if (local_ele < minimi && local_ele > 0.f)
									minimi = local_ele;
								d_E[idx] += local_ele;
							}
					}
						if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u)
							axACOSEM += (local_ele * d_OSEM[local_ind]);
#endif
#ifdef ATOMIC
						atomicAdd(&d_Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
						atomicAdd(&d_Summ[local_ind], local_ele);
#endif
#endif
					}
				}
			}
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
				const int dim_n = z_n * Nxy + y_n * (int)(Nx)+x_n;
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