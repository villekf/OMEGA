/**************************************************************************
* Special functions for the 3D orthogonal distance-based ray tracer.
*
* Copyright (C) 2020 Ville-Veikko Wettenhovi
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) ad_N1 later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT Ad_N1 WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <https://www.gnu.org/licenses/>.
***************************************************************************/
#pragma once
#include "general_orth_cuda_functions.cuh"

// Denominator (forward projection), multi-GPU version, orthogonal distance based ray tracer
__device__ void orth_distance_multi_psf(int tempi, const unsigned int d_N0, const unsigned int d_N4, float* temp, const unsigned int d_N2, int tempj, int tempk,
	const float local_sino, const float local_ele, const float* local_psf, const float* d_OSEM, const unsigned int d_Nxy,
	const bool no_norm, const bool RHS, const bool SUMMA, const unsigned int im_dim, CAST* Summ, const unsigned int d_N1, const unsigned int d_N3, const char
	start, int ju, int ku, unsigned char xyz, float* d_store_elements,
#ifdef MBSREM
	unsigned int* d_store_indices, unsigned int* ind, float* ax, const RecMethodsOpenCL MethodListOpenCL, const float d_h, float* d_E, const unsigned int idx,
	CAST* d_co, CAST* d_aco, float* minimi, const unsigned char MBSREM_prepass, float* axCOSEM, const unsigned int d_alku
#else
	unsigned int* d_store_indices, unsigned int* ind, CAST* d_rhs_OSEM, const unsigned char* MethodList, const float d_h, float* ax
#endif
) {
	bool breikki1 = false;
	bool breikki2 = false;
	if (start != 0) {
		if (start < 0) {
			ju *= -1;
			ku *= -1;
		}
		if (start == 21 || start == -21) {
			if (ju > 0) {
				if (xyz == 1u) {
					for (int uu = 0; uu >= -PSF_LIMIT; uu--) {
						int x = tempi + uu;
						if (x < 0 || x >= d_N0)
							break;
						for (int yy = 0; yy >= -PSF_LIMIT; yy--) {
							int y = tempj + yy;
							if (y < 0 || y >= d_N1)
								break;
							computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
								&breikki1, x, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, uu, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
								d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
								d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
							computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
								&breikki2, x, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, d_N4, uu, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
								d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
								d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
						}
					}
					for (int uu = 0; uu <= PSF_LIMIT; uu++) {
						int x = tempi + uu;
						if (x < 0 || x >= d_N0)
							break;
						for (int yy = 0; yy >= -PSF_LIMIT; yy--) {
							int y = tempj + yy;
							if (y < 0 || y >= d_N1)
								break;
							computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
								&breikki1, x, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, uu, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
								d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
								d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
							computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
								&breikki2, x, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, d_N4, uu, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
								d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
								d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
						}
					}
				}
				else {
					for (int yy = 0; yy >= -PSF_LIMIT; yy--) {
						int y = tempj + yy;
						if (y < 0 || y >= d_N1)
							break;
						for (int uu = 0; uu >= -PSF_LIMIT; uu--) {
							int x = tempi + uu;
							if (x < 0 || x >= d_N0)
								break;
							computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
								&breikki1, x, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, uu, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
								d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
								d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
							computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
								&breikki2, x, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, d_N4, uu, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
								d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
								d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
						}
						for (int uu = 0; uu <= PSF_LIMIT; uu++) {
							int x = tempi + uu;
							if (x < 0 || x >= d_N0)
								break;
							computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
								&breikki1, x, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, uu, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
								d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
								d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
							computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
								&breikki2, x, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, d_N4, uu, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
								d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
								d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
						}
					}
				}
			}
			else {
				if (xyz == 1u) {
					for (int uu = 0; uu >= -PSF_LIMIT; uu--) {
						int x = tempi + uu;
						if (x < 0 || x >= d_N0)
							break;
						for (int yy = 0; yy <= PSF_LIMIT; yy++) {
							int y = tempj + yy;
							if (y < 0 || y >= d_N1)
								break;
							computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
								&breikki1, x, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, uu, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
								d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
								d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
							computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
								&breikki2, x, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, d_N4, uu, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
								d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
								d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
						}
					}
					for (int uu = 0; uu <= PSF_LIMIT; uu++) {
						int x = tempi + uu;
						if (x < 0 || x >= d_N0)
							break;
						for (int yy = 0; yy <= PSF_LIMIT; yy++) {
							int y = tempj + yy;
							if (y < 0 || y >= d_N1)
								break;
							computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
								&breikki1, x, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, uu, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
								d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
								d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
							computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
								&breikki2, x, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, d_N4, uu, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
								d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
								d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
						}
					}
				}
				else {
					for (int yy = 0; yy <= PSF_LIMIT; yy++) {
						int y = tempj + yy;
						if (y < 0 || y >= d_N1)
							break;
						for (int uu = 0; uu >= -PSF_LIMIT; uu--) {
							int x = tempi + uu;
							if (x < 0 || x >= d_N0)
								break;
							computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
								&breikki1, x, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, uu, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
								d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
								d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
							computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
								&breikki2, x, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, d_N4, uu, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
								d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
								d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
						}
						for (int uu = 0; uu <= PSF_LIMIT; uu++) {
							int x = tempi + uu;
							if (x < 0 || x >= d_N0)
								break;
							computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
								&breikki1, x, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, uu, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
								d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
								d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
							computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
								&breikki2, x, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, d_N4, uu, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
								d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
								d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
						}
					}
				}
			}
		}
		else if (start == 22) {
			for (int uu = 0; uu >= -PSF_LIMIT; uu--) {
				int x = tempi + uu;
				if (x < 0 || x >= d_N0)
					break;
				computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
					&breikki1, x, tempj, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, uu, 0, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
					d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
					d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
				computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
					&breikki2, x, tempj, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, d_N4, uu, 0, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
					d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
					d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
			}
			for (int uu = 0; uu <= PSF_LIMIT; uu++) {
				int x = tempi + uu;
				if (x < 0 || x >= d_N0)
					break;
				computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
					&breikki1, x, tempj, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, uu, 0, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
					d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
					d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
				computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
					&breikki2, x, tempj, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, d_N4, uu, 0, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
					d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
					d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
			}
		}
		else if (start == 10) {
			if (xyz == 1u) {
				for (int uu = 0; uu >= -PSF_LIMIT; uu--) {
					int x = tempi + uu;
					if (x < 0 || x >= d_N0)
						break;
					computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
						&breikki1, x, tempj, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, uu, 0, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
						d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
						d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
						&breikki2, x, tempj, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, d_N4, uu, 0, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
						d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
						d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
				}
				for (int uu = 0; uu <= PSF_LIMIT; uu++) {
					int x = tempi + uu;
					if (x < 0 || x >= d_N0)
						break;
					computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
						&breikki1, x, tempj, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, uu, 0, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
						d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
						d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
						&breikki2, x, tempj, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, d_N4, uu, 0, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
						d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
						d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
				}
			}
			else {
				for (int yy = 0; yy >= -PSF_LIMIT; yy--) {
					int y = tempj + yy;
					if (y < 0 || y >= d_N1)
						break;
					computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
						&breikki1, tempi, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, 0, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
						d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
						d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
						&breikki2, tempi, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, d_N4, 0, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
						d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
						d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
				}
				for (int yy = 0; yy <= PSF_LIMIT; yy++) {
					int y = tempj + yy;
					if (y < 0 || y >= d_N1)
						break;
					computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
						&breikki1, tempi, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, 0, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
						d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
						d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
						&breikki2, tempi, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, d_N4, 0, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
						d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
						d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
				}
			}
		}
		else if (start == 11 || start == -11) {
			if (ju > 0) {
				if (xyz == 1u) {
					for (int uu = 0; uu >= -PSF_LIMIT; uu--) {
						int x = tempi + uu;
						if (x < 0 || x >= d_N0)
							break;
						computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
							&breikki1, x, tempj, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, uu, 0, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					}
					for (int uu = 0; uu <= PSF_LIMIT; uu++) {
						int x = tempi + uu;
						if (x < 0 || x >= d_N0)
							break;
						computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
							&breikki2, x, tempj, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, uu, 0, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					}
				}
				else {
					for (int yy = 0; yy >= -PSF_LIMIT; yy--) {
						int y = tempj + yy;
						if (y < 0 || y >= d_N1)
							break;
						computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
							&breikki1, tempi, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, 0, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
						computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
							&breikki2, tempi, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, d_N4, 0, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					}
				}
			}
			else {
				if (xyz == 1u) {
					for (int uu = 0; uu >= -PSF_LIMIT; uu--) {
						int x = tempi + uu;
						if (x < 0 || x >= d_N0)
							break;
						computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
							&breikki1, x, tempj, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, uu, 0, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					}
					for (int uu = 0; uu <= PSF_LIMIT; uu++) {
						int x = tempi + uu;
						if (x < 0 || x >= d_N0)
							break;
						computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
							&breikki2, x, tempj, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, d_N4, uu, 0, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					}
				}
				else {
					for (int yy = 0; yy <= PSF_LIMIT; yy++) {
						int y = tempj + yy;
						if (y < 0 || y >= d_N1)
							break;
						computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
							&breikki1, tempi, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, 0, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
						computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
							&breikki2, tempi, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, d_N4, 0, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					}
				}
			}
		}
		else if (start == 1 || start == -1) {
			if (ju > 0) {
				if (xyz == 1u) {
					for (int uu = 0; uu >= -PSF_LIMIT; uu--) {
						int x = tempi + uu;
						if (x < 0 || x >= d_N0)
							break;
						computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
							&breikki1, x, tempj, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, uu, 0, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					}
					for (int uu = 0; uu <= PSF_LIMIT; uu++) {
						int x = tempi + uu;
						if (x < 0 || x >= d_N0)
							break;
						computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
							&breikki2, x, tempj, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, uu, 0, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					}
				}
				else {
					for (int yy = 0; yy >= -PSF_LIMIT; yy--) {
						int y = tempj + yy;
						if (y < 0 || y >= d_N1)
							break;
						computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
							&breikki1, tempi, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, 0, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
						computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
							&breikki2, tempi, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, d_N4, 0, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					}
				}
			}
			else {
				if (xyz == 1u) {
					for (int uu = 0; uu >= -PSF_LIMIT; uu--) {
						int x = tempi + uu;
						if (x < 0 || x >= d_N0)
							break;
						computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
							&breikki1, x, tempj, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, d_N4, uu, 0, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					}
					for (int uu = 0; uu <= PSF_LIMIT; uu++) {
						int x = tempi + uu;
						if (x < 0 || x >= d_N0)
							break;
						computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
							&breikki2, x, tempj, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, d_N4, uu, 0, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					}
				}
				else {
					for (int yy = 0; yy <= PSF_LIMIT; yy++) {
						int y = tempj + yy;
						if (y < 0 || y >= d_N1)
							break;
						computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
							&breikki1, tempi, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, 0, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
						computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
							&breikki2, tempi, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, d_N4, 0, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					}
				}
			}
			if (ku > 0 && tempk > 0) {
				for (int zz = 0; zz >= -PSF_LIMIT; zz--) {
					int z = tempk + zz;
					if (z < 0 || z >= d_N4)
						break;
					for (int uu = 0; uu >= -PSF_LIMIT; uu--) {
						int x = tempi + uu;
						if (x < 0 || x >= d_N0)
							break;
						computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempk, tempj, local_sino, ax,
							&breikki1, x, z, no_norm, RHS, SUMMA, im_dim, d_N0, d_N4, d_N1, d_N2, d_Nxy, d_N3, uu, zz, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
						computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempk, tempj, local_sino, ax,
							&breikki2, x, z, no_norm, RHS, SUMMA, im_dim, d_N0, d_N4, d_N1, d_N2, d_Nxy, d_N3, d_N4, uu, zz, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					}
					for (int uu = 0; uu <= PSF_LIMIT; uu++) {
						int x = tempi + uu;
						if (x < 0 || x >= d_N0)
							break;
						computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempk, tempj, local_sino, ax,
							&breikki1, x, z, no_norm, RHS, SUMMA, im_dim, d_N0, d_N4, d_N1, d_N2, d_Nxy, d_N3, uu, zz, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
						computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempk, tempj, local_sino, ax,
							&breikki2, x, z, no_norm, RHS, SUMMA, im_dim, d_N0, d_N4, d_N1, d_N2, d_Nxy, d_N3, d_N4, uu, zz, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					}
				}
			}
			else {
				for (int zz = 0; zz <= PSF_LIMIT; zz++) {
					int z = tempk + zz;
					if (z < 0 || z >= d_N4)
						break;
					for (int uu = 0; uu >= -PSF_LIMIT; uu--) {
						int x = tempi + uu;
						if (x < 0 || x >= d_N0)
							break;
						computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempk, tempj, local_sino, ax,
							&breikki1, x, z, no_norm, RHS, SUMMA, im_dim, d_N0, d_N4, d_N1, d_N2, d_Nxy, d_N3, uu, zz, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
						computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempk, tempj, local_sino, ax,
							&breikki2, x, z, no_norm, RHS, SUMMA, im_dim, d_N0, d_N4, d_N1, d_N2, d_Nxy, d_N3, d_N4, uu, zz, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					}
					for (int uu = 0; uu <= PSF_LIMIT; uu++) {
						int x = tempi + uu;
						if (x < 0 || x >= d_N0)
							break;
						computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempk, tempj, local_sino, ax,
							&breikki1, x, z, no_norm, RHS, SUMMA, im_dim, d_N0, d_N4, d_N1, d_N2, d_Nxy, d_N3, uu, zz, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
						computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempk, tempj, local_sino, ax,
							&breikki2, x, z, no_norm, RHS, SUMMA, im_dim, d_N0, d_N4, d_N1, d_N2, d_Nxy, d_N3, d_N4, uu, zz, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					}
				}
			}
		}
		else if (start == 2 || start == -2) {
			if (xyz == 1u) {
				for (int uu = 0; uu >= -PSF_LIMIT; uu--) {
					int x = tempi + uu;
					if (x < 0 || x >= d_N0)
						break;
					computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
						&breikki1, x, tempj, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, uu, 0, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
						d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
						d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
						&breikki2, x, tempj, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, d_N4, uu, 0, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
						d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
						d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
				}
				for (int uu = 0; uu <= PSF_LIMIT; uu++) {
					int x = tempi + uu;
					if (x < 0 || x >= d_N0)
						break;
					computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
						&breikki1, x, tempj, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, uu, 0, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
						d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
						d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
						&breikki2, x, tempj, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, d_N4, uu, 0, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
						d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
						d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
				}
			}
			else {
				for (int yy = 0; yy >= -PSF_LIMIT; yy--) {
					int y = tempj + yy;
					if (y < 0 || y >= d_N1)
						break;
					computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
						&breikki1, tempi, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, 0, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
						d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
						d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
						&breikki2, tempi, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, d_N4, 0, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
						d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
						d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
				}
				for (int yy = 0; yy <= PSF_LIMIT; yy++) {
					int y = tempj + yy;
					if (y < 0 || y >= d_N1)
						break;
					computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
						&breikki1, tempi, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, 0, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
						d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
						d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempj, tempk, local_sino, ax,
						&breikki2, tempi, y, no_norm, RHS, SUMMA, im_dim, d_N0, d_N1, d_N4, d_N2, d_N3, d_Nxy, d_N4, 0, yy, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
						d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
						d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
				}
			}
			if (ku > 0) {
				for (int zz = 0; zz >= -PSF_LIMIT; zz--) {
					int z = tempk + zz;
					if (z < 0 || z >= d_N4)
						break;
					for (int uu = 0; uu >= -PSF_LIMIT; uu--) {
						int x = tempi + uu;
						if (x < 0 || x >= d_N0)
							break;
						computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempk, tempj, local_sino, ax,
							&breikki1, x, z, no_norm, RHS, SUMMA, im_dim, d_N0, d_N4, d_N1, d_N2, d_Nxy, d_N3, uu, zz, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
						computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempk, tempj, local_sino, ax,
							&breikki2, x, z, no_norm, RHS, SUMMA, im_dim, d_N0, d_N4, d_N1, d_N2, d_Nxy, d_N3, d_N4, uu, zz, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					}
					for (int uu = 0; uu <= PSF_LIMIT; uu++) {
						int x = tempi + uu;
						computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempk, tempj, local_sino, ax,
							&breikki1, x, z, no_norm, RHS, SUMMA, im_dim, d_N0, d_N4, d_N1, d_N2, d_Nxy, d_N3, uu, zz, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
						computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempk, tempj, local_sino, ax,
							&breikki2, x, z, no_norm, RHS, SUMMA, im_dim, d_N0, d_N4, d_N1, d_N2, d_Nxy, d_N3, d_N4, uu, zz, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					}
				}
			}
			else {
				for (int zz = 0; zz <= PSF_LIMIT; zz++) {
					int z = tempk + zz;
					if (z < 0 || z >= d_N4)
						break;
					for (int uu = 0; uu >= -PSF_LIMIT; uu--) {
						int x = tempi + uu;
						if (x < 0 || x >= d_N0)
							break;
						computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempk, tempj, local_sino, ax,
							&breikki1, x, z, no_norm, RHS, SUMMA, im_dim, d_N0, d_N4, d_N1, d_N2, d_Nxy, d_N3, uu, zz, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
						computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempk, tempj, local_sino, ax,
							&breikki2, x, z, no_norm, RHS, SUMMA, im_dim, d_N0, d_N4, d_N1, d_N2, d_Nxy, d_N3, d_N4, uu, zz, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					}
					for (int uu = 0; uu <= PSF_LIMIT; uu++) {
						int x = tempi + uu;
						if (x < 0 || x >= d_N0)
							break;
						computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempi, tempk, tempj, local_sino, ax,
							&breikki1, x, z, no_norm, RHS, SUMMA, im_dim, d_N0, d_N4, d_N1, d_N2, d_Nxy, d_N3, uu, zz, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
						computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempi, tempk, tempj, local_sino, ax,
							&breikki2, x, z, no_norm, RHS, SUMMA, im_dim, d_N0, d_N4, d_N1, d_N2, d_Nxy, d_N3, d_N4, uu, zz, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
							d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
							d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
					}
				}
			}
		}
	}
	else {
		computeOrthVoxelDecreasingPSF(local_ele, local_psf, temp, tempk, tempj, tempi, local_sino, ax,
			&breikki1, tempk, tempj, no_norm, RHS, SUMMA, im_dim, d_N4, d_N1, d_N0, d_Nxy, d_N3, d_N2, 0, 0, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
			d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
			d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
		computeOrthVoxelIncreasingPSF(local_ele, local_psf, temp, tempk, tempj, tempi, local_sino, ax,
			&breikki1, tempk, tempj, no_norm, RHS, SUMMA, im_dim, d_N4, d_N1, d_N0, d_Nxy, d_N3, d_N2, d_N4, 0, 0, Summ, d_OSEM, d_store_elements,
#ifdef MBSREM
			d_store_indices, ind, d_E, MethodListOpenCL, d_h, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
			d_store_indices, ind, d_rhs_OSEM, MethodList, d_h);
#endif
	}
}


__device__ void orth_distance_perpendicular_multi_3D(const float* center1, const float center2, const float* z_center,
	float* temp, float* ax, const float d_b, const float d, const float d_d1,
	const unsigned int d_N1, const unsigned int d_N2, const unsigned int z_loop, const float* d_atten, const float d_norm, const float local_sino, const unsigned int d_N,
	const unsigned int d_NN, const float* d_OSEM, const float xs, const float ys, const float zs, const float xl, const float yl, const float zl,
	const float crystal_size_z, const unsigned int d_N1x, const unsigned int d_N4, const unsigned char no_norm, CAST* Summ,
	const bool FPbool, const bool RHS, const float global_factor, const float bmin, const float bmax, const float Vmax, const float* V,
#ifdef MBSREM
	const RecMethodsOpenCL MethodListOpenCL, const unsigned int d_alku, float* axCOSEM, const float d_h,
	float* d_E, CAST* d_co, CAST* d_aco, float* minimi, const unsigned char MBSREM_prepass,
	const float* d_sc_ra, float* d_Amin, float* d_ACOSEM_lhs, const unsigned int idx
#else
	CAST* d_rhs_OSEM, const unsigned int im_dim, const unsigned char* MethodList, const float d_h
#endif
) {
	//const unsigned int zz = z_loop * d_N2 * d_N1;
	const unsigned int apu = perpendicular_start(d_b, d, d_d1, d_N1);
	float jelppi = 0.f;
	for (int zz = (int)(z_loop); zz >= 0; zz--) {
		for (int uu = (int)(apu); uu >= 0; uu--) {
#ifdef VOL
			float local_ele = compute_element_vol_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele >= bmax)
				break;
			if (local_ele < bmin)
				local_ele = Vmax;
			else
				local_ele = V[__float2uint_rn((local_ele - bmin) * CC)];
#else
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
#endif
			unsigned int local_ind = uu * d_N + zz * d_N1x;
			if (FPbool) {
				*temp += (local_ele * d_N2);
#if defined(ATN) || defined(MBSREM)
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
#ifdef ATN
					if (zz == (int)(z_loop) && uu == (int)(apu))
						jelppi += (d_d1 * -d_atten[local_ind]);
#endif
#ifdef MBSREM
					if (local_sino > 0.f && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0) {
						*axCOSEM += (local_ele * d_OSEM[local_ind]);
					}
#else
					if (local_sino > 0.f) {
						denominator(local_ele, ax, local_ind, im_dim, d_OSEM);
					}
#endif
					local_ind += d_NN;
				}
#endif
			}
			else if (RHS) {
				local_ele *= *temp;
#ifdef MBSREM
				if (d_alku == 0 && (MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
					if (local_ele < *minimi && local_ele > 0.f)
						* minimi = local_ele;
					d_E[idx] += local_ele;
				}
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (d_alku == 0) {
						if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
#ifdef ATOMIC
							atomicAdd(&d_co[local_ind], __float2ull_rn(*axCOSEM * (d_OSEM[local_ind] * local_ele) * TH));
#else
							atomicAdd(&d_co[local_ind], *axCOSEM * (d_OSEM[local_ind] * local_ele));
#endif
						if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
#ifdef ATOMIC
							atomicAdd(&d_aco[local_ind], __float2ull_rn(*axCOSEM * (powf(d_OSEM[local_ind], d_h) * local_ele) * TH));
#else
							atomicAdd(&d_aco[local_ind], *axCOSEM * (powf(d_OSEM[local_ind], d_h) * local_ele));
#endif
						if (MBSREM_prepass == 1)
#ifdef ATOMIC
							atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
							atomicAdd(&Summ[local_ind], local_ele);
#endif
					}
					else
						*ax += (local_ele * d_OSEM[local_ind]);
					local_ind += d_NN;
				}
#else
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (no_norm == 0u)
#ifdef ATOMIC
						atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
						atomicAdd(&Summ[local_ind], local_ele);
#endif
					rhs(MethodList, local_ele, ax, local_ind, im_dim, d_rhs_OSEM, d_h, d_OSEM);
					local_ind += d_NN;
				}
#endif
			}
			else {
				local_ele *= *temp;
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
#ifdef ATOMIC
					atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
					atomicAdd(&Summ[local_ind], local_ele);
#endif
					local_ind += d_NN;
				}
			}
		}
		for (unsigned int uu = apu + 1; uu < d_N1; uu++) {
#ifdef VOL
			float local_ele = compute_element_vol_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele >= bmax)
				break;
			if (local_ele < bmin)
				local_ele = Vmax;
			else
				local_ele = V[__float2uint_rn((local_ele - bmin) * CC)];
#else
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
#endif
			unsigned int local_ind = uu * d_N + zz * d_N1x;
			if (FPbool) {
				*temp += (local_ele * d_N2);
#ifdef MBSREM
				if (local_sino > 0.f && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0) {
					*axCOSEM += (local_ele * d_OSEM[local_ind]);
				}
#else
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (local_sino > 0.f) {
						denominator(local_ele, ax, local_ind, im_dim, d_OSEM);
					}
					local_ind += d_NN;
				}
#endif
			}
			else if (RHS) {
				local_ele *= *temp;
#ifdef MBSREM
				if (d_alku == 0 && (MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
					if (local_ele < *minimi && local_ele > 0.f)
						* minimi = local_ele;
					d_E[idx] += local_ele;
				}
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (d_alku == 0) {
						if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
#ifdef ATOMIC
							atomicAdd(&d_co[local_ind], __float2ull_rn(*axCOSEM * (d_OSEM[local_ind] * local_ele) * TH));
#else
							atomicAdd(&d_co[local_ind], *axCOSEM * (d_OSEM[local_ind] * local_ele));
#endif
						if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
#ifdef ATOMIC
							atomicAdd(&d_aco[local_ind], __float2ull_rn(*axCOSEM * (powf(d_OSEM[local_ind], d_h) * local_ele) * TH));
#else
							atomicAdd(&d_aco[local_ind], *axCOSEM * (powf(d_OSEM[local_ind], d_h) * local_ele));
#endif
						if (MBSREM_prepass == 1)
#ifdef ATOMIC
							atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
							atomicAdd(&Summ[local_ind], local_ele);
#endif
					}
					else
						*ax += (local_ele * d_OSEM[local_ind]);
					local_ind += d_NN;
				}
#else
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (no_norm == 0u)
#ifdef ATOMIC
						atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
						atomicAdd(&Summ[local_ind], local_ele);
#endif
					rhs(MethodList, local_ele, ax, local_ind, im_dim, d_rhs_OSEM, d_h, d_OSEM);
					local_ind += d_NN;
				}
#endif
			}
			else {
				local_ele *= *temp;
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
#ifdef ATOMIC
					atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
					atomicAdd(&Summ[local_ind], local_ele);
#endif
					local_ind += d_NN;
				}
			}
		}
	}
	for (unsigned int zz = z_loop + 1u; zz < d_N4; zz++) {
		for (int uu = (int)(apu); uu >= 0; uu--) {
#ifdef VOL
			float local_ele = compute_element_vol_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele >= bmax)
				break;
			if (local_ele < bmin)
				local_ele = Vmax;
			else
				local_ele = V[__float2uint_rn((local_ele - bmin) * CC)];
#else
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
#endif
			unsigned int local_ind = uu * d_N + zz * d_N1x;
			if (FPbool) {
				*temp += (local_ele * d_N2);
#ifdef MBSREM
				if (local_sino > 0.f && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0) {
					*axCOSEM += (local_ele * d_OSEM[local_ind]);
				}
#else
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (local_sino > 0.f) {
						denominator(local_ele, ax, local_ind, im_dim, d_OSEM);
					}
					local_ind += d_NN;
				}
#endif
			}
			else if (RHS) {
				local_ele *= *temp;
#ifdef MBSREM
				if (d_alku == 0 && (MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
					if (local_ele < *minimi && local_ele > 0.f)
						* minimi = local_ele;
					d_E[idx] += local_ele;
				}
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (d_alku == 0) {
						if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
#ifdef ATOMIC
							atomicAdd(&d_co[local_ind], __float2ull_rn(*axCOSEM * (d_OSEM[local_ind] * local_ele) * TH));
#else
							atomicAdd(&d_co[local_ind], *axCOSEM * (d_OSEM[local_ind] * local_ele));
#endif
						if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
#ifdef ATOMIC
							atomicAdd(&d_aco[local_ind], __float2ull_rn(*axCOSEM * (powf(d_OSEM[local_ind], d_h) * local_ele) * TH));
#else
							atomicAdd(&d_aco[local_ind], *axCOSEM * (powf(d_OSEM[local_ind], d_h) * local_ele));
#endif
						if (MBSREM_prepass == 1)
#ifdef ATOMIC
							atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
							atomicAdd(&Summ[local_ind], local_ele);
#endif
					}
					else
						*ax += (local_ele * d_OSEM[local_ind]);
					local_ind += d_NN;
				}
#else
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (no_norm == 0u)
#ifdef ATOMIC
						atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
						atomicAdd(&Summ[local_ind], local_ele);
#endif
					rhs(MethodList, local_ele, ax, local_ind, im_dim, d_rhs_OSEM, d_h, d_OSEM);
					local_ind += d_NN;
				}
#endif
			}
			else {
				local_ele *= *temp;
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
#ifdef ATOMIC
					atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
					atomicAdd(&Summ[local_ind], local_ele);
#endif
					local_ind += d_NN;
				}
			}
		}
		for (unsigned int uu = apu + 1; uu < d_N1; uu++) {
#ifdef VOL
			float local_ele = compute_element_vol_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele >= bmax)
				break;
			if (local_ele < bmin)
				local_ele = Vmax;
			else
				local_ele = V[__float2uint_rn((local_ele - bmin) * CC)];
#else
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
#endif
			unsigned int local_ind = uu * d_N + zz * d_N1x;
			if (FPbool) {
				*temp += (local_ele * d_N2);
#ifdef MBSREM
				if (local_sino > 0.f && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0) {
					*axCOSEM += (local_ele * d_OSEM[local_ind]);
				}
#else
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (local_sino > 0.f) {
						denominator(local_ele, ax, local_ind, im_dim, d_OSEM);
					}
					local_ind += d_NN;
				}
#endif
			}
			else if (RHS) {
				local_ele *= *temp;
#ifdef MBSREM
				if (d_alku == 0 && (MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
					if (local_ele < *minimi && local_ele > 0.f)
						* minimi = local_ele;
					d_E[idx] += local_ele;
				}
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (d_alku == 0) {
						if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
#ifdef ATOMIC
							atomicAdd(&d_co[local_ind], __float2ull_rn(*axCOSEM * (d_OSEM[local_ind] * local_ele) * TH));
#else
							atomicAdd(&d_co[local_ind], *axCOSEM * (d_OSEM[local_ind] * local_ele));
#endif
						if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
#ifdef ATOMIC
							atomicAdd(&d_aco[local_ind], __float2ull_rn(*axCOSEM * (powf(d_OSEM[local_ind], d_h) * local_ele) * TH));
#else
							atomicAdd(&d_aco[local_ind], *axCOSEM * (powf(d_OSEM[local_ind], d_h) * local_ele));
#endif
						if (MBSREM_prepass == 1)
#ifdef ATOMIC
							atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
							atomicAdd(&Summ[local_ind], local_ele);
#endif
					}
					else
						*ax += (local_ele * d_OSEM[local_ind]);
					local_ind += d_NN;
				}
#else
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (no_norm == 0u)
#ifdef ATOMIC
						atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
						atomicAdd(&Summ[local_ind], local_ele);
#endif
					rhs(MethodList, local_ele, ax, local_ind, im_dim, d_rhs_OSEM, d_h, d_OSEM);
					local_ind += d_NN;
				}
#endif
			}
			else {
				local_ele *= *temp;
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
#ifdef ATOMIC
					atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
					atomicAdd(&Summ[local_ind], local_ele);
#endif
					local_ind += d_NN;
				}
			}
		}
	}
#ifdef MBSREM
	if (!RHS) {
		*temp = 1.f / *temp;
#ifdef ATN
		* temp *= expf(jelppi);
#endif
#ifdef NORM
		* temp *= d_norm[idx];
#endif
		* temp *= global_factor;
	}
	else {
		if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1)
			d_Amin[idx] = *minimi;
		if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u) {
#ifdef RANDOMS
				* ax += d_sc_ra[idx];
#endif
			d_ACOSEM_lhs[idx] = *ax;
		}
	}
#else
	if (FPbool) {
		*temp = 1.f / *temp;
#ifdef ATN
		* temp *= expf(jelppi);
#endif
#ifdef NORM
		* temp *= d_norm;
#endif
		* temp *= global_factor;
	}
#endif
}
