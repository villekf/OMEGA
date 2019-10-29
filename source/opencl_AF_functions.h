/**************************************************************************
* Functions for implementation 2 kernel files. Uses f32 atomics.
*
* Copyright (C) 2019  Ville-Veikko Wettenhovi
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
#pragma once
#include "general_opencl_functions.h"
#include "opencl_AF_functions_nonatomic.h"

#define THR 0.001f



// Right-hand side of MLEM
inline void rhs_mlem(float local_ele, const __local float* ax, const uint local_ind, const uint lid, const uint n_rekos, 
	const uint d_N, __global float* d_rhs_MLEM) {
	uint yy = local_ind;
	for (uint kk = lid * n_rekos; kk < (lid + 1u) * n_rekos; kk++) {
		atomicAdd_g_f(&d_rhs_MLEM[yy], (local_ele * ax[kk]));
		yy += d_N;
	}
}

// Denominator (forward projection), orthogonal distance based ray tracer
inline void orth_distance_mlem(const int tempi, const uint Nx, const float y_diff, const float x_diff, float y_center, __constant float* x_center, 
	const float kerroin, const float length_, float* temp, const uint temp_ijk, const float local_sino, __local float* ax, const uint d_Ny, const uint d_N, 
	const uchar no_norm, const uint lid, const uint n_rekos, const uint im_dim, const __global float* d_MLEM, __global float* d_rhs_MLEM, __global float* Summ,
	const bool RHS, const bool SUMMA) {

	const float diff = kerroin - x_diff * y_center;
	for (int uu = tempi; uu >= 0; uu--) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		if (RHS) {
			local_ele *= *temp;
			rhs_mlem(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_rhs_MLEM);
				
			if (no_norm == 0u)
				atomicAdd_g_f(&Summ[local_ind], local_ele);
		}
		else if (SUMMA) {
			local_ele *= *temp;
			atomicAdd_g_f(&Summ[local_ind], local_ele);
		}
		else {
			*temp += local_ele;
			if (local_sino > 0.f) {
				denominator(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_MLEM);
			}
		}
	}
	for (int uu = tempi + 1; uu < convert_int(Nx); uu++) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		if (RHS) {
			local_ele *= *temp;
			rhs_mlem(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_rhs_MLEM);
				
			if (no_norm == 0u)
				atomicAdd_g_f(&Summ[local_ind], local_ele);
		}
		else if (SUMMA) {
			local_ele *= *temp;
			atomicAdd_g_f(&Summ[local_ind], local_ele);
		}
		else {
			*temp += local_ele;
			if (local_sino > 0.f) {
				denominator(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_MLEM);
			}
		}
	}
}

// Denominator (forward projection), multi-GPU version, orthogonal distance based ray tracer
inline void orth_distance_mlem_3D(const int tempi, const uint Nx, const uint Nz, const float y_diff, const float x_diff, const float z_diff,
	const float y_center, __constant float* x_center, __constant float* z_center, const float crystal_size_z, float* temp, const uint temp_ijk, 
	const float local_sino, __local float* ax, const uint d_Ny, const uint d_N, const int tempk, const uint Nxy, const float xs, const float ys, 
	const float zs, const int dec, const uchar no_norm, const uint lid, const uint n_rekos, const uint im_dim, const __global float* d_MLEM, 
	__global float* d_rhs_MLEM, __global float* Summ, const bool RHS, const bool SUMMA) {

	bool loppu1 = true;
	bool loppu2 = true;
	uchar loppu3 = 0u;
	uchar loppu4 = 0u;
	bool pass1 = true;
	bool pass = false;
	int alku = tempk;
	int alku2 = tempk;
	int alku3 = tempk;
	int alku4 = tempk + 1;
	for (int uu = tempi; uu >= 0; uu--) {
		//for (int zz = tempk; zz >= 0; zz--) {
		for (int zz = alku; zz >= 0; zz--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
			if (local_ele <= THR) {
				if (loppu1 && zz >= alku3 - dec && !(pass && alku3 == tempk)) {
					continue;
				}
				else {
					if (loppu1 && loppu3 == 1u && loppu4 == 0u)
						pass1 = false;
					break;
				}
			}
			const uint local_ind = compute_ind_orth_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
			if (RHS) {
				local_ele *= *temp;
				rhs_mlem(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_rhs_MLEM);
					
				if (no_norm == 0u)
					atomicAdd_g_f(&Summ[local_ind], local_ele);
			}
			else if (SUMMA) {
				local_ele *= *temp;
				atomicAdd_g_f(&Summ[local_ind], local_ele);
			}
			else {
				*temp += local_ele;
				if (local_sino > 0.f) {
					denominator(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_MLEM);
				}
			}
			if (loppu1) {
				alku = zz;
				loppu1 = false;
				loppu3 = 1u;
				if (zz == tempk || loppu4 == 1u)
					pass = true;
			}
			alku3 = zz;
		}
		if (pass) {
			for (int zz = alku2 + 1; zz < Nz; zz++) {
				//for (uint zz = convert_uint_sat(tempk) + 1u; zz < Nz; zz++) {
				float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
				if (local_ele <= THR) {
					if (loppu2 && zz <= alku4 + dec) {
						continue;
					}
					else {
						if (loppu2 && loppu4 == 1u) {
							loppu4 = 0u;
							pass = false;
						}
						break;
					}
				}
				const uint local_ind = compute_ind_orth_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
				if (RHS) {
					local_ele *= *temp;
					rhs_mlem(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_rhs_MLEM);
						
					if (no_norm == 0u)
						atomicAdd_g_f(&Summ[local_ind], local_ele);
				}
				else if (SUMMA) {
					local_ele *= *temp;
					atomicAdd_g_f(&Summ[local_ind], local_ele);
				}
				else {
					*temp += local_ele;
					if (local_sino > 0.f) {
						denominator(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_MLEM);
					}
				}
				if (loppu2) {
					alku2 = zz - 1;
					loppu2 = false;
					loppu4 = 1u;
					pass1 = true;
				}
				alku4 = zz;
			}
		}
		loppu1 = true;
		loppu2 = true;
		if (!pass1) {
			break;
		}
	}
	loppu3 = 0u;
	pass1 = true;
	pass = false;
	loppu4 = 0u;
	alku2 = tempk;
	alku = tempk;
	alku3 = tempk;
	alku4 = tempk + 1;
	for (uint uu = convert_uint_sat(tempi) + 1u; uu < Nx; uu++) {
		for (int zz = alku; zz >= 0; zz--) {
			//for (int zz = tempk; zz >= 0; zz--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
			if (local_ele <= THR) {
				if (loppu1 && zz >= alku3 - dec && !(pass && alku3 == tempk)) {
					continue;
				}
				else {
					if (loppu1 && loppu3 == 1u && loppu4 == 0u)
						pass1 = false;
					break;
				}
			}
			const uint local_ind = compute_ind_orth_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
			if (RHS) {
				local_ele *= *temp;
				rhs_mlem(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_rhs_MLEM);
					
				if (no_norm == 0u)
					atomicAdd_g_f(&Summ[local_ind], local_ele);
			}
			else if (SUMMA) {
				local_ele *= *temp;
				atomicAdd_g_f(&Summ[local_ind], local_ele);
			}
			else {
				*temp += local_ele;
				if (local_sino > 0.f) {
					denominator(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_MLEM);
				}
			}
			if (loppu1) {
				alku = zz;
				loppu1 = false;
				loppu3 = 1u;
				if (zz == tempk || loppu4 == 1u)
					pass = true;
			}
			alku3 = zz;
		}
		if (pass) {
			for (int zz = alku2 + 1; zz < Nz; zz++) {
				//for (uint zz = convert_uint_sat(tempk) + 1u; zz < Nz; zz++) {
				float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
				if (local_ele <= THR) {
					if (loppu2 && zz <= alku4 + dec) {
						continue;
					}
					else {
						if (loppu2 && loppu4 == 1u) {
							loppu4 = 0u;
							pass = false;
						}
						break;
					}
				}
				const uint local_ind = compute_ind_orth_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
				if (RHS) {
					local_ele *= *temp;
					rhs_mlem(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_rhs_MLEM);
						
					if (no_norm == 0u)
						atomicAdd_g_f(&Summ[local_ind], local_ele);
				}
				else if (SUMMA) {
					local_ele *= *temp;
					atomicAdd_g_f(&Summ[local_ind], local_ele);
				}
				else {
					*temp += local_ele;
					if (local_sino > 0.f) {
						denominator(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_MLEM);
					}
				}
				if (loppu2) {
					alku2 = zz - 1;
					loppu2 = false;
					loppu4 = 1u;
					pass1 = true;
				}
				alku4 = zz;
			}
		}
		loppu1 = true;
		loppu2 = true;
		if (!pass1) {
			break;
		}
	}
}

// Denominator (forward projection), orthogonal distance based ray tracer
inline void orth_distance_perpendicular_mlem(const float diff1, __constant float* center1, const float kerroin,
	const float length_, float* temp, const uint d_attenuation_correction, const uint d_normalization, __local float* ax,
	const float d_b, const float d, const float d_d1, const uint d_N1, const uint d_N2, const uint z_loop, 
	const __global float* d_atten, const __global float* d_norm, const uint idx, const float local_sino, const uint d_N, const uint d_NN, const uchar no_norm,
	const uint lid, const uint n_rekos, const uint im_dim, const __global float* d_MLEM, __global float* d_rhs_MLEM, __global float* Summ, const bool FP, const bool RHS) {

	const uint zz = z_loop * d_N2 * d_N1;
	const int apu = perpendicular_start(d_b, d, d_d1, d_N1);
	float jelppi = 0.f;
	for (int uu = apu; uu >= 0; uu--) {
		float local_ele = compute_element_orth(kerroin, diff1, center1[uu], length_);
		if (local_ele <= THR)
			break;
		uint local_ind = convert_uint_sat(uu) * d_N + zz;
		if (FP) {
			*temp += (local_ele * d_N2);
			for (uint kk = 0u; kk < d_N2; kk++) {
				if (d_attenuation_correction && uu == apu)
					jelppi += (d_d1 * -d_atten[local_ind]);
				if (local_sino > 0.f) {
					denominator(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_MLEM);
						
				}
				local_ind += d_NN;
			}
		}
		else if (RHS) {
			local_ele *= *temp;
			for (uint kk = 0u; kk < d_N2; kk++) {
				if (no_norm == 0u)
					atomicAdd_g_f(&Summ[local_ind], local_ele);
				rhs_mlem(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_rhs_MLEM);
					
				local_ind += d_NN;
			}
		}
		else {
			local_ele *= *temp;
			for (uint kk = 0u; kk < d_N2; kk++) {
				atomicAdd_g_f(&Summ[local_ind], local_ele);
				local_ind += d_NN;
			}
		}
	}
	for (uint uu = convert_uint_sat(apu) + 1u; uu < d_N1; uu++) {
		float local_ele = compute_element_orth(kerroin, diff1, center1[uu], length_);
		if (local_ele <= THR)
			break;
		uint local_ind = uu * d_N + zz;
		if (FP) {
			*temp += (local_ele * d_N2);
			for (uint kk = 0u; kk < d_N2; kk++) {
				if (local_sino > 0.f) {
					denominator(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_MLEM);
						
				}
				local_ind += d_NN;
			}
		}
		else if (RHS) {
			local_ele *= *temp;
			for (uint kk = 0u; kk < d_N2; kk++) {
				if (no_norm == 0u)
					atomicAdd_g_f(&Summ[local_ind], local_ele);
				rhs_mlem(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_rhs_MLEM);
					
			}
		}
		else {
			local_ele *= *temp;
			for (uint kk = 0u; kk < d_N2; kk++) {
				atomicAdd_g_f(&Summ[local_ind], local_ele);
				local_ind += d_NN;
			}
		}
	}
	if (FP) {
		*temp = 1.f / *temp;
		if (d_attenuation_correction)
			* temp *= jelppi;
		if (d_normalization == 1u)
			* temp *= d_norm[idx];
	}
}


inline void orth_distance_perpendicular_mlem_3D(__constant float* center1, const float center2, __constant float* z_center, float* temp, 
	const uint d_attenuation_correction, const uint d_normalization, __local float* ax, const float d_b, const float d, 
	const float d_d1, const uint d_N1, const uint d_N2, const uint z_loop, const __global float* d_atten, const __global float* d_norm, const uint idx, 
	const float local_sino, const uint d_N, const uint d_NN, const float xs, const float ys, const float zs, const float xl, const float yl, const float zl, 
	const float crystal_size_z, const uint Nyx, const uint Nz, const uchar no_norm, const uint lid, const uint n_rekos, const uint im_dim, 
	const __global float* d_MLEM, __global float* d_rhs_MLEM, __global float* Summ, const bool FP, const bool RHS) {

	//const uint zz = z_loop * d_N2 * d_N1;
	const uint apu = perpendicular_start(d_b, d, d_d1, d_N1);
	float jelppi = 0.f;
	for (int zz = convert_int_sat(z_loop); zz >= 0; zz--) {
		for (int uu = convert_int_sat(apu); uu >= 0; uu--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			if (FP) {
				*temp += (local_ele * d_N2);
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (d_attenuation_correction && zz == convert_int_sat(z_loop) && uu == convert_int_sat(apu))
						jelppi += (d_d1 * -d_atten[local_ind]);
					if (local_sino > 0.f) {
						denominator(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_MLEM);
							
					}
					local_ind += d_NN;
				}
			}
			else if (RHS) {
				local_ele *= *temp;
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (no_norm == 0u)
						atomicAdd_g_f(&Summ[local_ind], local_ele);
					rhs_mlem(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_rhs_MLEM);
						
					local_ind += d_NN;
				}
			}
			else {
				local_ele *= *temp;
				for (uint kk = 0u; kk < d_N2; kk++) {
					atomicAdd_g_f(&Summ[local_ind], local_ele);
					local_ind += d_NN;
				}
			}
		}
		for (uint uu = apu + 1; uu < d_N1; uu++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			if (FP) {
				*temp += (local_ele * d_N2);
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (local_sino > 0.f) {
						denominator(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_MLEM);
							
					}
					local_ind += d_NN;
				}
			}
			else if (RHS) {
				local_ele *= *temp;
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (no_norm == 0u)
						atomicAdd_g_f(&Summ[local_ind], local_ele);
					rhs_mlem(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_rhs_MLEM);
						
					local_ind += d_NN;
				}
			}
			else {
				local_ele *= *temp;
				for (uint kk = 0u; kk < d_N2; kk++) {
					atomicAdd_g_f(&Summ[local_ind], local_ele);
					local_ind += d_NN;
				}
			}
		}
	}
	for (uint zz = z_loop + 1u; zz < Nz; zz++) {
		for (int uu = convert_int_sat(apu); uu >= 0; uu--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			if (FP) {
				*temp += (local_ele * d_N2);
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (local_sino > 0.f) {
						denominator(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_MLEM);
							
					}
					local_ind += d_NN;
				}
			}
			else if (RHS) {
				local_ele *= *temp;
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (no_norm == 0u)
						atomicAdd_g_f(&Summ[local_ind], local_ele);
					rhs_mlem(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_rhs_MLEM);
						
					local_ind += d_NN;
				}
			}
			else {
				local_ele *= *temp;
				for (uint kk = 0u; kk < d_N2; kk++) {
					atomicAdd_g_f(&Summ[local_ind], local_ele);
					local_ind += d_NN;
				}
			}
		}
		for (uint uu = apu + 1; uu < d_N1; uu++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			if (FP) {
				*temp += (local_ele * d_N2);
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (local_sino > 0.f) {
						denominator(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_MLEM);
							
					}
					local_ind += d_NN;
				}
			}
			else if (RHS) {
				local_ele *= *temp;
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (no_norm == 0u)
						atomicAdd_g_f(&Summ[local_ind], local_ele);
					rhs_mlem(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_rhs_MLEM);
						
					local_ind += d_NN;
				}
			}
			else {
				local_ele *= *temp;
				for (uint kk = 0u; kk < d_N2; kk++) {
					atomicAdd_g_f(&Summ[local_ind], local_ele);
					local_ind += d_NN;
				}
			}
		}
	}
	if (FP) {
		*temp = 1. / *temp;
		if (d_attenuation_correction)
			* temp *= jelppi;
		if (d_normalization == 1u)
			* temp *= d_norm[idx];
	}
}

// Right-hand side
inline void rhs(__constant uchar* MethodList, float local_ele, const __local float* ax, const uint local_ind, const uint lid, const uint n_rekos,
	const uint d_N, __global float* d_rhs_OSEM, const float d_h, const __global float* d_ACOSEM) {
	uint yy = local_ind;
	uint oo = 0u;
	for (uint kk = lid * n_rekos; kk < (lid + 1u) * n_rekos; kk++) {
		if (MethodList[oo] < 2u)
			atomicAdd_g_f(&d_rhs_OSEM[yy], (local_ele * ax[kk]));
		else if (MethodList[oo] == 2u) // COSEM
			atomicAdd_g_f(&d_rhs_OSEM[yy], (ax[kk] * (local_ele * d_ACOSEM[yy])));
		else if (MethodList[oo] == 3u) //ACOSEM
			atomicAdd_g_f(&d_rhs_OSEM[yy], ax[kk] * (local_ele * native_powr(d_ACOSEM[yy], d_h)));
		yy += d_N;
		oo++;
	}
}

// RHS, orthogonal distance based ray tracer
inline void orth_distance_rhs(const int tempi, const uint Nx, const float y_diff, const float x_diff, const float y_center, __constant float* x_center, 
	const float kerroin, const float length_, const float temp, const uint temp_ijk, __constant uchar* MethodList, const __local float* ax, 
	const uint d_Ny, const uint d_N, const uchar no_norm, const uint lid, const uint n_rekos, const uint im_dim, __global float* d_rhs_OSEM, const float d_h, 
	const __global float* d_ACOSEM, __global float* Summ) {

	const float diff = kerroin - x_diff * y_center;
	for (int uu = tempi; uu >= 0; uu--) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		local_ele *= temp;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		if (no_norm == 0u)
			atomicAdd_g_f(&Summ[local_ind], local_ele);
		rhs(MethodList, local_ele, ax, local_ind, lid, n_rekos, im_dim, d_rhs_OSEM, d_h, d_ACOSEM);
	}
	for (int uu = tempi + 1; uu < Nx; uu++) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		local_ele *= temp;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		rhs(MethodList, local_ele, ax, local_ind, lid, n_rekos, im_dim, d_rhs_OSEM, d_h, d_ACOSEM);
		if (no_norm == 0u)
			atomicAdd_g_f(&Summ[local_ind], local_ele);
	}
}

// RHS, orthogonal distance based ray tracer
inline void orth_distance_summ(const int tempi, const uint Nx, const float y_diff, const float x_diff, const float y_center, __constant float* x_center, 
	const float kerroin, const float length_, const float temp, const uint temp_ijk, __global float* Summ, const uint d_Ny, const uint d_N) {

	const float diff = kerroin - x_diff * y_center;
	for (int uu = tempi; uu >= 0; uu--) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		local_ele *= temp;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		atomicAdd_g_f(&Summ[local_ind], local_ele);
	}
	for (int uu = tempi + 1; uu < Nx; uu++) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		local_ele *= temp;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		atomicAdd_g_f(&Summ[local_ind], local_ele);
	}
}

// RHS, orthogonal distance based ray tracer
inline void orth_distance_rhs_3D(const int tempi, const uint Nx, const uint Nz, const float y_diff, const float x_diff, const float z_diff,
	const float y_center, __constant float* x_center, __constant float* z_center, const float crystal_size_z, const float temp, const uint temp_ijk, 
	__constant uchar* MethodList, __local float* ax, const uint d_Ny, const uint d_N, const uint Nxy, const float xs, const float ys, const float zs, 
	const int tempk, const int dec, const uchar no_norm, const uint lid, const uint n_rekos, const uint im_dim,
	__global float* d_rhs_OSEM, const float d_h, const __global float* d_ACOSEM, __global float* Summ) {

	bool loppu1 = true;
	bool loppu2 = true;
	uchar loppu3 = 0u;
	uchar loppu4 = 0u;
	bool pass1 = true;
	bool pass = false;
	int alku = tempk;
	int alku2 = tempk;
	int alku3 = tempk;
	int alku4 = tempk + 1;
	for (int uu = tempi; uu >= 0; uu--) {
		//for (int zz = tempk; zz >= 0; zz--) {
		for (int zz = alku; zz >= 0; zz--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
			if (local_ele <= THR) {
				if (loppu1 && zz >= alku3 - dec && !(pass && alku3 == tempk)) {
					continue;
				}
				else {
					if (loppu1 && loppu3 == 1u && loppu4 == 0u)
						pass1 = false;
					break;
				}
			}
			const uint local_ind = compute_ind_orth_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
			local_ele *= temp;
			if (no_norm == 0u)
				atomicAdd_g_f(&Summ[local_ind], local_ele);
			rhs(MethodList, local_ele, ax, local_ind, lid, n_rekos, im_dim, d_rhs_OSEM, d_h, d_ACOSEM);
			if (loppu1) {
				alku = zz;
				loppu1 = false;
				loppu3 = 1u;
				if (zz == tempk || loppu4 == 1u)
					pass = true;
			}
			alku3 = zz;
		}
		if (pass) {
			for (int zz = alku2 + 1; zz < Nz; zz++) {
				//for (uint zz = convert_uint_sat(tempk) + 1u; zz < Nz; zz++) {
				float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
				if (local_ele <= THR) {
					if (loppu2 && zz <= alku4 + dec) {
						continue;
					}
					else {
						if (loppu2 && loppu4 == 1u) {
							loppu4 = 0u;
							pass = false;
						}
						break;
					}
				}
				const uint local_ind = compute_ind_orth_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
				local_ele *= temp;
				if (no_norm == 0u)
					atomicAdd_g_f(&Summ[local_ind], local_ele);
				rhs(MethodList, local_ele, ax, local_ind, lid, n_rekos, im_dim, d_rhs_OSEM, d_h, d_ACOSEM);
				if (loppu2) {
					alku2 = zz - 1;
					loppu2 = false;
					loppu4 = 1u;
					pass1 = true;
				}
				alku4 = zz;
			}
		}
		loppu1 = true;
		loppu2 = true;
		if (!pass1) {
			break;
		}
	}
	loppu3 = 0u;
	pass1 = true;
	pass = false;
	loppu4 = 0u;
	alku2 = tempk;
	alku = tempk;
	alku3 = tempk;
	alku4 = tempk + 1;
	for (uint uu = convert_uint_sat(tempi) + 1u; uu < Nx; uu++) {
		for (int zz = alku; zz >= 0; zz--) {
			//for (int zz = tempk; zz >= 0; zz--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
			if (local_ele <= THR) {
				if (loppu1 && zz >= alku3 - dec && !(pass && alku3 == tempk)) {
					continue;
				}
				else {
					if (loppu1 && loppu3 == 1u && loppu4 == 0u)
						pass1 = false;
					break;
				}
			}
			const uint local_ind = compute_ind_orth_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
			local_ele *= temp;
			if (no_norm == 0u)
				atomicAdd_g_f(&Summ[local_ind], local_ele);
			rhs(MethodList, local_ele, ax, local_ind, lid, n_rekos, im_dim, d_rhs_OSEM, d_h, d_ACOSEM);
			if (loppu1) {
				alku = zz;
				loppu1 = false;
				loppu3 = 1u;
				if (zz == tempk || loppu4 == 1u)
					pass = true;
			}
			alku3 = zz;
		}
		if (pass) {
			for (int zz = alku2 + 1; zz < Nz; zz++) {
				//for (uint zz = convert_uint_sat(tempk) + 1u; zz < Nz; zz++) {
				float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
				if (local_ele <= THR) {
					if (loppu2 && zz <= alku4 + dec) {
						continue;
					}
					else {
						if (loppu2 && loppu4 == 1u) {
							loppu4 = 0u;
							pass = false;
						}
						break;
					}
				}
				const uint local_ind = compute_ind_orth_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
				local_ele *= temp;
				if (no_norm == 0u)
					atomicAdd_g_f(&Summ[local_ind], local_ele);
				rhs(MethodList, local_ele, ax, local_ind, lid, n_rekos, im_dim, d_rhs_OSEM, d_h, d_ACOSEM);
				if (loppu2) {
					alku2 = zz - 1;
					loppu2 = false;
					loppu4 = 1u;
					pass1 = true;
				}
				alku4 = zz;
			}
		}
		loppu1 = true;
		loppu2 = true;
		if (!pass1) {
			break;
		}
	}
}

// RHS, orthogonal distance based ray tracer
inline void orth_distance_summ_3D(const int tempi, const uint Nx, const uint Nz, const float y_diff, const float x_diff, const float z_diff,
	const float y_center, __constant float* x_center, __constant float* z_center, const float crystal_size_z, const float temp, const uint temp_ijk,
	const uint d_Ny, const uint d_N, const uint Nxy, const float xs, const float ys, const float zs, const int tempk,
	const int dec, __global float* Summ) {

	bool loppu1 = true;
	bool loppu2 = true;
	uchar loppu3 = 0u;
	uchar loppu4 = 0u;
	bool pass1 = true;
	bool pass = false;
	int alku = tempk;
	int alku2 = tempk;
	int alku3 = tempk;
	int alku4 = tempk + 1;
	for (int uu = tempi; uu >= 0; uu--) {
		//for (int zz = tempk; zz >= 0; zz--) {
		for (int zz = alku; zz >= 0; zz--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
			if (local_ele <= THR) {
				if (loppu1 && zz >= alku3 - dec && !(pass && alku3 == tempk)) {
					continue;
				}
				else {
					if (loppu1 && loppu3 == 1u && loppu4 == 0u)
						pass1 = false;
					break;
				}
			}
			const uint local_ind = compute_ind_orth_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
			local_ele *= temp;
			atomicAdd_g_f(&Summ[local_ind], local_ele);
			if (loppu1) {
				alku = zz;
				loppu1 = false;
				loppu3 = 1u;
				if (zz == tempk || loppu4 == 1u)
					pass = true;
			}
			alku3 = zz;
		}
		if (pass) {
			for (int zz = alku2 + 1; zz < Nz; zz++) {
				//for (uint zz = convert_uint_sat(tempk) + 1u; zz < Nz; zz++) {
				float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
				if (local_ele <= THR) {
					if (loppu2 && zz <= alku4 + dec) {
						continue;
					}
					else {
						if (loppu2 && loppu4 == 1u) {
							loppu4 = 0u;
							pass = false;
						}
						break;
					}
				}
				const uint local_ind = compute_ind_orth_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
				local_ele *= temp;
				atomicAdd_g_f(&Summ[local_ind], local_ele);
				if (loppu2) {
					alku2 = zz - 1;
					loppu2 = false;
					loppu4 = 1u;
					pass1 = true;
				}
				alku4 = zz;
			}
		}
		loppu1 = true;
		loppu2 = true;
		if (!pass1) {
			break;
		}
	}
	loppu3 = 0u;
	pass1 = true;
	pass = false;
	loppu4 = 0u;
	alku2 = tempk;
	alku = tempk;
	alku3 = tempk;
	alku4 = tempk + 1;
	for (uint uu = convert_uint_sat(tempi) + 1u; uu < Nx; uu++) {
		for (int zz = alku; zz >= 0; zz--) {
			//for (int zz = tempk; zz >= 0; zz--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
			if (local_ele <= THR) {
				if (loppu1 && zz >= alku3 - dec && !(pass && alku3 == tempk)) {
					continue;
				}
				else {
					if (loppu1 && loppu3 == 1u && loppu4 == 0u)
						pass1 = false;
					break;
				}
			}
			const uint local_ind = compute_ind_orth_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
			local_ele *= temp;
			atomicAdd_g_f(&Summ[local_ind], local_ele);
			if (loppu1) {
				alku = zz;
				loppu1 = false;
				loppu3 = 1u;
				if (zz == tempk || loppu4 == 1u)
					pass = true;
			}
			alku3 = zz;
		}
		if (pass) {
			for (int zz = alku2 + 1; zz < Nz; zz++) {
				//for (uint zz = convert_uint_sat(tempk) + 1u; zz < Nz; zz++) {
				float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
				if (local_ele <= THR) {
					if (loppu2 && zz <= alku4 + dec) {
						continue;
					}
					else {
						if (loppu2 && loppu4 == 1u) {
							loppu4 = 0u;
							pass = false;
						}
						break;
					}
				}
				const uint local_ind = compute_ind_orth_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
				local_ele *= temp;
				atomicAdd_g_f(&Summ[local_ind], local_ele);
				if (loppu2) {
					alku2 = zz - 1;
					loppu2 = false;
					loppu4 = 1u;
					pass1 = true;
				}
				alku4 = zz;
			}
		}
		loppu1 = true;
		loppu2 = true;
		if (!pass1) {
			break;
		}
	}
}


// RHS, orthogonal distance based ray tracer
inline void orth_distance_rhs_perpendicular(const float diff1, __constant float* center1, const float kerroin, const float length_, const float temp, 
	__constant uchar* MethodList, __local float* ax,	const float d_b, const float d, const float d_d1, const uint d_N1, const uint d_N2, 
	const uint z_loop, const uint d_N, const uint d_NN, const uchar no_norm, const uint lid, const uint n_rekos, const uint im_dim, 
	__global float* d_rhs_OSEM,	const float d_h, const __global float* d_ACOSEM, __global float* Summ) {

	const uint zz = z_loop * d_N2 * d_N1;
	const int apu = perpendicular_start(d_b, d, d_d1, d_N1);
	for (int uu = apu; uu >= 0; uu--) {
		float local_ele = compute_element_orth(kerroin, diff1, center1[uu], length_);
		if (local_ele <= THR)
			break;
		local_ele *= temp;
		uint local_ind = convert_uint_sat(uu) * d_N + zz;
		for (uint kk = 0u; kk < d_N2; kk++) {
			//atomicAdd_g_f(&d_rhs_OSEM[local_ind], 1.f);
			if (no_norm == 0u)
				atomicAdd_g_f(&Summ[local_ind], local_ele);
			rhs(MethodList, local_ele, ax, local_ind, lid, n_rekos, im_dim, d_rhs_OSEM, d_h, d_ACOSEM);
			local_ind += d_NN;
		}
	}
	for (uint uu = convert_uint_sat(apu) + 1u; uu < d_N1; uu++) {
		float local_ele = compute_element_orth(kerroin, diff1, center1[uu], length_);
		if (local_ele <= THR)
			break;
		local_ele *= temp;
		uint local_ind = uu * d_N + zz;
		for (uint kk = 0u; kk < d_N2; kk++) {
			//atomicAdd_g_f(&d_rhs_OSEM[local_ind], 1.f);
			rhs(MethodList, local_ele, ax, local_ind, lid, n_rekos, im_dim, d_rhs_OSEM, d_h, d_ACOSEM);
			if (no_norm == 0u)
				atomicAdd_g_f(&Summ[local_ind], local_ele);
			local_ind += d_NN;
		}
	}
}

// Compute the normalizer, orthogonal distance based ray tracer
inline void orth_distance_summ_perpendicular(const float diff1, __constant float* center1, const float kerroin, const float length_, const float temp, 
	const float d_b, const float d, const float d_d1, const uint d_N1, const uint d_N2, const uint z_loop, const uint d_N, const uint d_NN,
	__global float* Summ) {

	const uint zz = z_loop * d_N2 * d_N1;
	const int apu = perpendicular_start(d_b, d, d_d1, d_N1);
	for (int uu = apu; uu >= 0; uu--) {
		float local_ele = compute_element_orth(kerroin, diff1, center1[uu], length_);
		if (local_ele <= THR)
			break;
		local_ele *= temp;
		uint local_ind = convert_uint_sat(uu) * d_N + zz;
		for (uint kk = 0u; kk < d_N2; kk++) {
			atomicAdd_g_f(&Summ[local_ind], local_ele);
			local_ind += d_NN;
		}
	}
	for (uint uu = convert_uint_sat(apu) + 1u; uu < d_N1; uu++) {
		float local_ele = compute_element_orth(kerroin, diff1, center1[uu], length_);
		if (local_ele <= THR)
			break;
		local_ele *= temp;
		uint local_ind = uu * d_N + zz;
		for (uint kk = 0u; kk < d_N2; kk++) {
			atomicAdd_g_f(&Summ[local_ind], local_ele);
			local_ind += d_NN;
		}
	}
}

inline void orth_distance_rhs_perpendicular_3D(__constant float* center1, const float center2, __constant float* z_center, const float temp, 
	__constant uchar* MethodList, __local float* ax,const float d_b, const float d, const float d_d1, const uint d_N1, const uint d_N2, const uint z_loop,
	const uint d_N, const uint d_NN, const uchar no_norm, const float xs, const float ys, const float zs, const float xl, const float yl, const float zl, 
	const float crystal_size_z, const uint Nyx, const uint Nz, const uint lid, const uint n_rekos, const uint im_dim, __global float* d_rhs_OSEM,
	const float d_h, const __global float* d_ACOSEM, __global float* Summ) {

	const uint apu = perpendicular_start(d_b, d, d_d1, d_N1);
	float jelppi = 0.f;
	for (int zz = convert_int_sat(z_loop); zz >= 0; zz--) {
		for (int uu = convert_int_sat(apu); uu >= 0; uu--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			local_ele *= temp;
			for (uint kk = 0u; kk < d_N2; kk++) {
				if (no_norm == 0u)
					atomicAdd_g_f(&Summ[local_ind], local_ele);
				rhs(MethodList, local_ele, ax, local_ind, lid, n_rekos, im_dim, d_rhs_OSEM, d_h, d_ACOSEM);
				local_ind += d_NN;
			}
		}
		for (uint uu = apu + 1; uu < d_N1; uu++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			local_ele *= temp;
			for (uint kk = 0u; kk < d_N2; kk++) {
				if (no_norm == 0u)
					atomicAdd_g_f(&Summ[local_ind], local_ele);
				rhs(MethodList, local_ele, ax, local_ind, lid, n_rekos, im_dim, d_rhs_OSEM, d_h, d_ACOSEM);
				local_ind += d_NN;
			}
		}
	}
	for (uint zz = z_loop + 1u; zz < Nz; zz++) {
		for (int uu = convert_int_sat(apu); uu >= 0; uu--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			local_ele *= temp;
			for (uint kk = 0u; kk < d_N2; kk++) {
				if (no_norm == 0u)
					atomicAdd_g_f(&Summ[local_ind], local_ele);
				rhs(MethodList, local_ele, ax, local_ind, lid, n_rekos, im_dim, d_rhs_OSEM, d_h, d_ACOSEM);
				local_ind += d_NN;
			}
		}
		for (uint uu = apu + 1; uu < d_N1; uu++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			local_ele *= temp;
			for (uint kk = 0u; kk < d_N2; kk++) {
				if (no_norm == 0u)
					atomicAdd_g_f(&Summ[local_ind], local_ele);
				rhs(MethodList, local_ele, ax, local_ind, lid, n_rekos, im_dim, d_rhs_OSEM, d_h, d_ACOSEM);
				local_ind += d_NN;
			}
		}
	}
}


inline void orth_distance_summ_perpendicular_3D(__constant float* center1, const float center2, __constant float* z_center, const float temp, 
	const float d_b, const float d, const float d_d1, const uint d_N1, const uint d_N2, const uint z_loop, const uint d_N, const uint d_NN, 
	const float xs, const float ys, const float zs, const float xl, const float yl, const float zl, const float crystal_size_z, const uint Nyx,
	const uint Nz, __global float* Summ) {

	const uint apu = perpendicular_start(d_b, d, d_d1, d_N1);
	float jelppi = 0.f;
	for (int zz = convert_int_sat(z_loop); zz >= 0; zz--) {
		for (int uu = convert_int_sat(apu); uu >= 0; uu--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			local_ele *= temp;
			for (uint kk = 0u; kk < d_N2; kk++) {
				atomicAdd_g_f(&Summ[local_ind], local_ele);
				local_ind += d_NN;
			}
		}
		for (uint uu = apu + 1; uu < d_N1; uu++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			local_ele *= temp;
			for (uint kk = 0u; kk < d_N2; kk++) {
				atomicAdd_g_f(&Summ[local_ind], local_ele);
				local_ind += d_NN;
			}
		}
	}
	for (uint zz = z_loop + 1u; zz < Nz; zz++) {
		for (int uu = convert_int_sat(apu); uu >= 0; uu--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			local_ele *= temp;
			for (uint kk = 0u; kk < d_N2; kk++) {
				atomicAdd_g_f(&Summ[local_ind], local_ele);
				local_ind += d_NN;
			}
		}
		for (uint uu = apu + 1; uu < d_N1; uu++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			local_ele *= temp;
			for (uint kk = 0u; kk < d_N2; kk++) {
				atomicAdd_g_f(&Summ[local_ind], local_ele);
				local_ind += d_NN;
			}
		}
	}
}

inline void orth_distance_perpendicular_cosem_3D(__constant float* center1, const float center2, __constant float* z_center, float* temp, 
	const uint d_attenuation_correction, const uint d_normalization, const float d_b, const float d, const float d_d1, const uint d_N1, const uint d_N2, 
	const uint z_loop, const __global float* d_atten, const __global float* d_norm, const uint idx, const uint d_N, const float local_sino,
	const uint d_NN, const float xs, const float ys, const float zs, const float xl, const float yl, const float zl, const float crystal_size_z, const uint Nyx, 
	const uint Nz, const RecMethodsOpenCL MethodListOpenCL, const uint d_alku, float* axCOSEM, const __global float* d_COSEM, const float d_h,
	__global float* d_E, __global float* d_co, __global float* d_aco, float* minimi, const uchar MBSREM_prepass, __global float* d_Summ, float* axACOSEM,
	const __global float* d_sc_ra, const uint d_randoms, __global float* d_Amin, __global float* d_ACOSEM_lhs, const bool RHS) {

	const uint apu = perpendicular_start(d_b, d, d_d1, d_N1);
	float jelppi = 0.f;
	for (int zz = convert_int_sat(z_loop); zz >= 0; zz--) {
		for (int uu = convert_int_sat(apu); uu >= 0; uu--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			if (RHS) {
				local_ele *= *temp;
				if (d_alku == 0 && (MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1) {
					if (local_ele < *minimi && local_ele != 0.f)
						* minimi = local_ele;
					d_E[idx] += local_ele;
				}
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (d_alku == 0) {
						if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
							atomicAdd_g_f(&d_co[local_ind], (*axCOSEM * local_ele * d_COSEM[local_ind]));
						if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
							atomicAdd_g_f(&d_aco[local_ind], *axCOSEM * (native_powr(d_COSEM[local_ind], d_h) * local_ele));
						if (MBSREM_prepass == 1)
							atomicAdd_g_f(&d_Summ[local_ind], local_ele);
					}
					else
						*axACOSEM += (local_ele * d_COSEM[local_ind]);
					local_ind += d_NN;
				}
			}
			else {
				*temp += local_ele;
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (d_attenuation_correction && zz == convert_int_sat(z_loop) && uu == convert_int_sat(apu))
						jelppi += (d_d1 * -d_atten[local_ind]);
					if (local_sino > 0.f && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0) {
						*axCOSEM += (local_ele * d_COSEM[local_ind]);
					}
					local_ind += d_NN;
				}
			}
		}
		for (uint uu = apu + 1; uu < d_N1; uu++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			if (RHS) {
				local_ele *= *temp;
				if (d_alku == 0 && (MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1) {
					if (local_ele < *minimi && local_ele != 0.f)
						* minimi = local_ele;
					d_E[idx] += local_ele;
				}
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (d_alku == 0) {
						if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
							atomicAdd_g_f(&d_co[local_ind], (*axCOSEM * local_ele * d_COSEM[local_ind]));
						if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
							atomicAdd_g_f(&d_aco[local_ind], *axCOSEM * (native_powr(d_COSEM[local_ind], d_h) * local_ele));
						if (MBSREM_prepass == 1)
							atomicAdd_g_f(&d_Summ[local_ind], local_ele);
					}
					else
						*axACOSEM += (local_ele * d_COSEM[local_ind]);
					local_ind += d_NN;
				}
			}
			else {
				*temp += local_ele;
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (local_sino > 0.f && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0) {
						*axCOSEM += (local_ele * d_COSEM[local_ind]);
					}
					local_ind += d_NN;
				}
			}
		}
	}
	for (uint zz = z_loop + 1u; zz < Nz; zz++) {
		for (int uu = convert_int_sat(apu); uu >= 0; uu--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			if (RHS) {
				local_ele *= *temp;
				if (d_alku == 0 && (MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1) {
					if (local_ele < *minimi && local_ele != 0.f)
						* minimi = local_ele;
					d_E[idx] += local_ele;
				}
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (d_alku == 0) {
						if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
							atomicAdd_g_f(&d_co[local_ind], (*axCOSEM * local_ele * d_COSEM[local_ind]));
						if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
							atomicAdd_g_f(&d_aco[local_ind],* axCOSEM * (native_powr(d_COSEM[local_ind], d_h) * local_ele));
						if (MBSREM_prepass == 1)
							atomicAdd_g_f(&d_Summ[local_ind], local_ele);
					}
					else
						*axACOSEM += (local_ele * d_COSEM[local_ind]);
					local_ind += d_NN;
				}
			}
			else {
				*temp += local_ele;
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (local_sino > 0.f && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0) {
						*axCOSEM += (local_ele * d_COSEM[local_ind]);
					}
					local_ind += d_NN;
				}
			}
		}
		for (uint uu = apu + 1; uu < d_N1; uu++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			if (RHS) {
				local_ele *= *temp;
				if (d_alku == 0 && (MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1) {
					if (local_ele < *minimi && local_ele != 0.f)
						* minimi = local_ele;
					d_E[idx] += local_ele;
				}
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (d_alku == 0) {
						if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
							atomicAdd_g_f(&d_co[local_ind], (*axCOSEM * local_ele * d_COSEM[local_ind]));
						if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
							atomicAdd_g_f(&d_aco[local_ind], *axCOSEM * (native_powr(d_COSEM[local_ind], d_h) * local_ele));
						if (MBSREM_prepass == 1)
							atomicAdd_g_f(&d_Summ[local_ind], local_ele);
					}
					else
						*axACOSEM += (local_ele * d_COSEM[local_ind]);
					local_ind += d_NN;
				}
			}
			else {
				*temp += local_ele;
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (local_sino > 0.f && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0) {
						*axCOSEM += (local_ele * d_COSEM[local_ind]);
					}
					local_ind += d_NN;
				}
			}
		}
	}
	if (!RHS) {
		*temp = 1.f / *temp;
		if (d_attenuation_correction == 1u)
			* temp *= jelppi;
		if (d_normalization == 1u)
			* temp *= d_norm[idx];
	}
	else {
		if ((MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1)
			d_Amin[idx] = *minimi;
		if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u) {
			if (d_randoms == 1u)
				*axACOSEM += d_sc_ra[idx];
			d_ACOSEM_lhs[idx] = *axACOSEM;
		}
	}
}

inline void orth_distance_perpendicular_cosem(const float diff1, __constant float* center1, const float kerroin, const float length_, 
	float* temp, const uint d_attenuation_correction, const uint d_normalization, const float d_b, const float d, const float d_d1, const uint d_N1, 
	const uint d_N2, const uint z_loop, const __global float* d_atten, const __global float* d_norm, const uint idx, const uint d_N, const float local_sino,
	const uint d_NN, const RecMethodsOpenCL MethodListOpenCL, const uint d_alku, float* axCOSEM, const __global float* d_COSEM, const float d_h, 
	__global float* d_E, __global float* d_co, __global float* d_aco, float* minimi, const uchar MBSREM_prepass, __global float* d_Summ, float* axACOSEM, 
	const __global float* d_sc_ra, const uint d_randoms, __global float* d_Amin, __global float* d_ACOSEM_lhs, const bool RHS) {

	const uint zz = z_loop * d_N2 * d_N1;
	const int apu = perpendicular_start(d_b, d, d_d1, d_N1);
	float jelppi = 0.f;
	for (int uu = apu; uu >= 0; uu--) {
		float local_ele = compute_element_orth(kerroin, diff1, center1[uu], length_);
		if (local_ele <= THR)
			break;
		uint local_ind = convert_uint_sat(uu) * d_N + zz;
		if (RHS) {
			local_ele *= *temp;
			if (d_alku == 0 && (MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1) {
				if (local_ele < *minimi && local_ele != 0.f)
					* minimi = local_ele;
				d_E[idx] += local_ele;
			}
			for (uint kk = 0u; kk < d_N2; kk++) {
				if (d_alku == 0) {
					if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
						atomicAdd_g_f(&d_co[local_ind], (*axCOSEM * local_ele * d_COSEM[local_ind]));
					if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
						atomicAdd_g_f(&d_aco[local_ind], *axCOSEM * (native_powr(d_COSEM[local_ind], d_h) * local_ele));
					if (MBSREM_prepass == 1)
						atomicAdd_g_f(&d_Summ[local_ind], local_ele);
				}
				else
					*axACOSEM += (local_ele * d_COSEM[local_ind]);
				local_ind += d_NN;
			}
		}
		else {
			*temp += (local_ele * d_N2);
			for (uint kk = 0u; kk < d_N2; kk++) {
				if (d_attenuation_correction && uu == apu)
					jelppi += (d_d1 * -d_atten[local_ind]);
				if (local_sino > 0.f && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0) {
					*axCOSEM += (local_ele * d_COSEM[local_ind]);
				}
				local_ind += d_NN;
			}
		}
	}
	for (uint uu = convert_uint_sat(apu) + 1u; uu < d_N1; uu++) {
		float local_ele = compute_element_orth(kerroin, diff1, center1[uu], length_);
		if (local_ele <= THR)
			break;
		uint local_ind = uu * d_N + zz;
		if (RHS) {
			local_ele *= *temp;
			if (d_alku == 0 && (MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1) {
				if (local_ele < *minimi && local_ele != 0.f)
					* minimi = local_ele;
				d_E[idx] += local_ele;
			}
			for (uint kk = 0u; kk < d_N2; kk++) {
				if (d_alku == 0) {
					if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
						atomicAdd_g_f(&d_co[local_ind], (*axCOSEM * local_ele * d_COSEM[local_ind]));
					if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
						atomicAdd_g_f(&d_aco[local_ind], *axCOSEM * (native_powr(d_COSEM[local_ind], d_h) * local_ele));
					if (MBSREM_prepass == 1)
						atomicAdd_g_f(&d_Summ[local_ind], local_ele);
				}
				else
					*axACOSEM += (local_ele * d_COSEM[local_ind]);
				local_ind += d_NN;
			}
		}
		else {
			*temp += (local_ele * d_N2);
			for (uint kk = 0u; kk < d_N2; kk++) {
				if (local_sino > 0.f && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0) {
					*axCOSEM += (local_ele * d_COSEM[local_ind]);
				}
				local_ind += d_NN;
			}
		}
	}
	if (RHS) {
		if ((MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1)
			d_Amin[idx] = *minimi;
		if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u) {
			if (d_randoms == 1u)
				* axACOSEM += d_sc_ra[idx];
			d_ACOSEM_lhs[idx] = *axACOSEM;
		}
	}
	else {
		*temp = 1.f / *temp;
		if (d_attenuation_correction == 1u)
			* temp *= jelppi;
		if (d_normalization == 1u)
			* temp *= d_norm[idx];
	}
}

// Denominator (forward projection), orthogonal distance based ray tracer
inline void orth_distance_cosem(const int tempi, const uint Nx, const float y_diff, const float x_diff, const float y_center, __constant float* x_center, 
	const float kerroin, const float length_, float* temp, const uint temp_ijk, const RecMethodsOpenCL MethodListOpenCL, const uint d_alku,
	const float local_sino, const uint d_Ny, const uint d_N, float* axCOSEM, const __global float* d_COSEM, const float d_h, __global float* d_E, const uint idx, 
	__global float* d_co, __global float* d_aco, float* minimi, const uchar MBSREM_prepass, __global float* d_Summ, float* axACOSEM, const bool RHS) {

	const float diff = kerroin - x_diff * y_center;
	for (int uu = tempi; uu >= 0; uu--) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		if (RHS) {
			local_ele *= *temp;
			if (d_alku == 0 && (MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1) {
				if (local_ele < *minimi && local_ele != 0.f)
					* minimi = local_ele;
				d_E[idx] += local_ele;
			}
			if (d_alku == 0) {
				if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
					atomicAdd_g_f(&d_co[local_ind], (*axCOSEM * local_ele * d_COSEM[local_ind]));
				if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
					atomicAdd_g_f(&d_aco[local_ind], *axCOSEM * (native_powr(d_COSEM[local_ind], d_h) * local_ele));
				if (MBSREM_prepass == 1)
					atomicAdd_g_f(&d_Summ[local_ind], local_ele);
			}
			else
				*axACOSEM += (local_ele * d_COSEM[local_ind]);
		}
		else {
			*temp += local_ele;
			if (local_sino > 0.f && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0) {
				*axCOSEM += (local_ele * d_COSEM[local_ind]);
			}
		}
	}
	for (int uu = tempi + 1; uu < convert_int(Nx); uu++) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		if (RHS) {
			local_ele *= *temp;
			if (d_alku == 0 && (MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1) {
				if (local_ele < *minimi && local_ele != 0.f)
					* minimi = local_ele;
				d_E[idx] += local_ele;
			}
			if (d_alku == 0) {
				if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
					atomicAdd_g_f(&d_co[local_ind], (*axCOSEM * local_ele * d_COSEM[local_ind]));
				if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
					atomicAdd_g_f(&d_aco[local_ind], *axCOSEM * (native_powr(d_COSEM[local_ind], d_h) * local_ele));
				if (MBSREM_prepass == 1)
					atomicAdd_g_f(&d_Summ[local_ind], local_ele);
			}
			else
				*axACOSEM += (local_ele * d_COSEM[local_ind]);
		}
		else {
			*temp += local_ele;
			if (local_sino > 0.f && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0) {
				*axCOSEM += (local_ele * d_COSEM[local_ind]);
			}
		}
	}
}


// Denominator (forward projection), orthogonal distance based ray tracer
inline void orth_distance_cosem_3D(const int tempi, const uint Nx, const uint Nz, const float y_diff, const float x_diff, const float z_diff,
	const float y_center, __constant float* x_center, __constant float* z_center, const float crystal_size_z, float* temp, const uint temp_ijk, 
	const RecMethodsOpenCL MethodListOpenCL, const uint d_alku, const float local_sino, const uint d_Ny, const uint d_N, const int tempk, const uint Nxy, 
	const float xs, const float ys, const float zs, float* axCOSEM, const __global float* d_COSEM, const float d_h, __global float* d_E, const uint idx, 
	__global float* d_co, __global float* d_aco, float* minimi, const uchar MBSREM_prepass, __global float* d_Summ, float* axACOSEM, const int dec, const bool RHS) {

	bool loppu1 = true;
	bool loppu2 = true;
	uchar loppu3 = 0u;
	uchar loppu4 = 0u;
	bool pass1 = true;
	bool pass = false;
	int alku = tempk;
	int alku2 = tempk;
	int alku3 = tempk;
	int alku4 = tempk + 1;
	for (int uu = tempi; uu >= 0; uu--) {
		//for (int zz = tempk; zz >= 0; zz--) {
		for (int zz = alku; zz >= 0; zz--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
			if (local_ele <= THR) {
				if (loppu1 && zz >= alku3 - dec && !(pass && alku3 == tempk)) {
					continue;
				}
				else {
					if (loppu1 && loppu3 == 1u && loppu4 == 0u)
						pass1 = false;
					break;
				}
			}
			const uint local_ind = compute_ind_orth_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
			if (RHS) {
				local_ele *= *temp;
				if (d_alku == 0u && (MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1) {
					if (local_ele < *minimi && local_ele != 0.f)
						* minimi = local_ele;
					d_E[idx] += local_ele;
				}
				if (d_alku == 0u) {
					if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
						atomicAdd_g_f(&d_co[local_ind], (*axCOSEM * local_ele * d_COSEM[local_ind]));
					if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
						atomicAdd_g_f(&d_aco[local_ind], *axCOSEM * (native_powr(d_COSEM[local_ind], d_h) * local_ele));
					if (MBSREM_prepass == 1)
						atomicAdd_g_f(&d_Summ[local_ind], local_ele);
				}
				else
					*axACOSEM += (local_ele * d_COSEM[local_ind]);
			}
			else {
				*temp += local_ele;
				if (local_sino > 0.f && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0) {
					*axCOSEM += (local_ele * d_COSEM[local_ind]);
				}
			}
			if (loppu1) {
				alku = zz;
				loppu1 = false;
				loppu3 = 1u;
				if (zz == tempk || loppu4 == 1u)
					pass = true;
			}
			alku3 = zz;
		}
		if (pass) {
			for (int zz = alku2 + 1; zz < Nz; zz++) {
				//for (uint zz = convert_uint_sat(tempk) + 1u; zz < Nz; zz++) {
				float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
				if (local_ele <= THR) {
					if (loppu2 && zz <= alku4 + dec) {
						continue;
					}
					else {
						if (loppu2 && loppu4 == 1u) {
							loppu4 = 0u;
							pass = false;
						}
						break;
					}
				}
				const uint local_ind = compute_ind_orth_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
				if (RHS) {
					local_ele *= *temp;
					if (d_alku == 0u && (MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1) {
						if (local_ele < *minimi && local_ele != 0.f)
							* minimi = local_ele;
						d_E[idx] += local_ele;
					}
					if (d_alku == 0u) {
						if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
							atomicAdd_g_f(&d_co[local_ind], (*axCOSEM * local_ele * d_COSEM[local_ind]));
						if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
							atomicAdd_g_f(&d_aco[local_ind], *axCOSEM * (native_powr(d_COSEM[local_ind], d_h) * local_ele));
						if (MBSREM_prepass == 1)
							atomicAdd_g_f(&d_Summ[local_ind], local_ele);
					}
					else
						*axACOSEM += (local_ele * d_COSEM[local_ind]);
				}
				else {
					*temp += local_ele;
					if (local_sino > 0.f && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0) {
						*axCOSEM += (local_ele * d_COSEM[local_ind]);
					}
				}
				if (loppu2) {
					alku2 = zz - 1;
					loppu2 = false;
					loppu4 = 1u;
					pass1 = true;
				}
				alku4 = zz;
			}
		}
		loppu1 = true;
		loppu2 = true;
		if (!pass1) {
			break;
		}
	}
	loppu3 = 0u;
	pass1 = true;
	pass = false;
	loppu4 = 0u;
	alku2 = tempk;
	alku = tempk;
	alku3 = tempk;
	alku4 = tempk + 1;
	for (uint uu = convert_uint_sat(tempi) + 1u; uu < Nx; uu++) {
		for (int zz = alku; zz >= 0; zz--) {
			//for (int zz = tempk; zz >= 0; zz--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
			if (local_ele <= THR) {
				if (loppu1 && zz >= alku3 - dec && !(pass && alku3 == tempk)) {
					continue;
				}
				else {
					if (loppu1 && loppu3 == 1u && loppu4 == 0u)
						pass1 = false;
					break;
				}
			}
			const uint local_ind = compute_ind_orth_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
			if (RHS) {
				local_ele *= *temp;
				if (d_alku == 0u && (MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1) {
					if (local_ele < *minimi && local_ele != 0.f)
						* minimi = local_ele;
					d_E[idx] += local_ele;
				}
				if (d_alku == 0u) {
					if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
						atomicAdd_g_f(&d_co[local_ind], (*axCOSEM * local_ele * d_COSEM[local_ind]));
					if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
						atomicAdd_g_f(&d_aco[local_ind], *axCOSEM * (native_powr(d_COSEM[local_ind], d_h) * local_ele));
					if (MBSREM_prepass == 1)
						atomicAdd_g_f(&d_Summ[local_ind], local_ele);
				}
				else
					*axACOSEM += (local_ele * d_COSEM[local_ind]);
			}
			else {
				*temp += local_ele;
				if (local_sino > 0.f && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0) {
					*axCOSEM += (local_ele * d_COSEM[local_ind]);
				}
			}
			if (loppu1) {
				alku = zz;
				loppu1 = false;
				loppu3 = 1u;
				if (zz == tempk || loppu4 == 1u)
					pass = true;
			}
			alku3 = zz;
		}
		if (pass) {
			for (int zz = alku2 + 1; zz < Nz; zz++) {
				//for (uint zz = convert_uint_sat(tempk) + 1u; zz < Nz; zz++) {
				float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
				if (local_ele <= THR) {
					if (loppu2 && zz <= alku4 + dec) {
						continue;
					}
					else {
						if (loppu2 && loppu4 == 1u) {
							loppu4 = 0u;
							pass = false;
						}
						break;
					}
				}
				const uint local_ind = compute_ind_orth_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
				if (RHS) {
					local_ele *= *temp;
					if (d_alku == 0u && (MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1) {
						if (local_ele < *minimi && local_ele != 0.f)
							* minimi = local_ele;
						d_E[idx] += local_ele;
					}
					if (d_alku == 0u) {
						if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
							atomicAdd_g_f(&d_co[local_ind], (*axCOSEM * local_ele * d_COSEM[local_ind]));
						if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
							atomicAdd_g_f(&d_aco[local_ind], *axCOSEM * (native_powr(d_COSEM[local_ind], d_h) * local_ele));
						if (MBSREM_prepass == 1)
							atomicAdd_g_f(&d_Summ[local_ind], local_ele);
					}
					else
						*axACOSEM += (local_ele * d_COSEM[local_ind]);
				}
				else {
					*temp += local_ele;
					if (local_sino > 0.f && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0) {
						*axCOSEM += (local_ele * d_COSEM[local_ind]);
					}
				}
				if (loppu2) {
					alku2 = zz - 1;
					loppu2 = false;
					loppu4 = 1u;
					pass1 = true;
				}
				alku4 = zz;
			}
		}
		loppu1 = true;
		loppu2 = true;
		if (!pass1) {
			break;
		}
	}
}