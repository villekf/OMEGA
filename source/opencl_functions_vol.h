/**************************************************************************
* Special functions for the volume-based ray tracer.
*
* Copyright (C) 2020  Ville-Veikko Wettenhovi
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

inline void vol_perpendicular_multi_3D(__constant float* center1, const float center2, __constant float* z_center,
	float* temp, const uint d_attenuation_correction, const uint d_normalization, float* ax, const float d_b, const float d, const float d_d1,
	const uint d_N1, const uint d_N2, const uint z_loop, const __global float* d_atten, const float d_norm, const float local_sino, const uint d_N, 
	const uint d_NN, const __global float* d_OSEM, const float xs, const float ys, const float zs, const float xl, const float yl, const float zl, 
	const float crystal_size_z, const uint Nyx, const uint Nz, const uchar no_norm, __global CAST* Summ, __global CAST* d_rhs_OSEM, 
	const bool FP_bool, const bool RHS, const float global_factor, const float bmin, const float bmax, const float Vmax, __constant float* V) {

	//const uint zz = z_loop * d_N2 * d_N1;
	const uint apu = perpendicular_start(d_b, d, d_d1, d_N1);
	float jelppi = 0.f;
	for (int zz = convert_int_sat(z_loop); zz >= 0; zz--) {
		for (int uu = convert_int_sat(apu); uu >= 0; uu--) {
			float local_ele = compute_element_vol_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele >= bmax)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			if (local_ele < bmin)
				local_ele = Vmax;
			else
				local_ele = V[convert_uint_rte((local_ele - bmin) * CC)];
			if (FP_bool) {
				*temp += (local_ele * d_N2);
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (d_attenuation_correction && zz == convert_int_sat(z_loop) && uu == convert_int_sat(apu))
						jelppi += (d_d1 * -d_atten[local_ind]);
					if (local_sino > 0.f) {
						denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
					}
					local_ind += d_NN;
				}
			}
			else if (RHS) {
				local_ele *= *temp;
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (no_norm == 0u)
#ifdef ATOMIC
						atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
#else
						atomicAdd_g_f(&Summ[local_ind], local_ele);
#endif
#ifdef ATOMIC
					atom_add(&d_rhs_OSEM[local_ind], convert_ulong_sat(local_ele * *ax * TH));
#else
					atomicAdd_g_f(&d_rhs_OSEM[local_ind], (local_ele * *ax));
#endif
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
			float local_ele = compute_element_vol_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele >= bmax)
				break;
			if (local_ele < bmin)
				local_ele = Vmax;
			else
				local_ele = V[convert_uint_rte((local_ele - bmin) * CC)];
			uint local_ind = uu * d_N + zz * Nyx;
			if (FP_bool) {
				*temp += (local_ele * d_N2);
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (local_sino > 0.f) {
						denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
					}
					local_ind += d_NN;
				}
			}
			else if (RHS) {
				local_ele *= *temp;
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (no_norm == 0u)
#ifdef ATOMIC
						atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
#else
						atomicAdd_g_f(&Summ[local_ind], local_ele);
#endif
#ifdef ATOMIC
					atom_add(&d_rhs_OSEM[local_ind], convert_ulong_sat(local_ele * *ax * TH));
#else
					atomicAdd_g_f(&d_rhs_OSEM[local_ind], (local_ele * *ax));
#endif
					local_ind += d_NN;
				}
			}
			else {
				local_ele *= *temp;
				for (uint kk = 0u; kk < d_N2; kk++) {
#ifdef ATOMIC
					atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
#else
					atomicAdd_g_f(&Summ[local_ind], local_ele);
#endif
					local_ind += d_NN;
				}
			}
		}
	}
	for (uint zz = z_loop + 1u; zz < Nz; zz++) {
		for (int uu = convert_int_sat(apu); uu >= 0; uu--) {
			float local_ele = compute_element_vol_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele >= bmax)
				break;
			if (local_ele < bmin)
				local_ele = Vmax;
			else
				local_ele = V[convert_uint_rte((local_ele - bmin) * CC)];
			uint local_ind = uu * d_N + zz * Nyx;
			if (FP_bool) {
				*temp += (local_ele * d_N2);
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (local_sino > 0.f) {
						denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
					}
					local_ind += d_NN;
				}
			}
			else if (RHS) {
				local_ele *= *temp;
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (no_norm == 0u)
#ifdef ATOMIC
						atom_add(&Summ[local_ind], convert_ulong_sat(local_ele* TH));
#else
						atomicAdd_g_f(&Summ[local_ind], local_ele);
#endif
#ifdef ATOMIC
					atom_add(&d_rhs_OSEM[local_ind], convert_ulong_sat(local_ele * *ax* TH));
#else
					atomicAdd_g_f(&d_rhs_OSEM[local_ind], (local_ele * *ax));
#endif
					local_ind += d_NN;
				}
			}
			else {
				local_ele *= *temp;
				for (uint kk = 0u; kk < d_N2; kk++) {
#ifdef ATOMIC
					atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
#else
					atomicAdd_g_f(&Summ[local_ind], local_ele);
#endif
					local_ind += d_NN;
				}
			}
		}
		for (uint uu = apu + 1; uu < d_N1; uu++) {
			float local_ele = compute_element_vol_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele >= bmax)
				break;
			if (local_ele < bmin)
				local_ele = Vmax;
			else
				local_ele = V[convert_uint_rte((local_ele - bmin) * CC)];
			uint local_ind = uu * d_N + zz * Nyx;
			if (FP_bool) {
				*temp += (local_ele * d_N2);
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (local_sino > 0.f) {
						denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
					}
					local_ind += d_NN;
				}
			}
			else if (RHS) {
				local_ele *= *temp;
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (no_norm == 0u)
#ifdef ATOMIC
						atom_add(&Summ[local_ind], convert_ulong_sat(local_ele* TH));
#else
						atomicAdd_g_f(&Summ[local_ind], local_ele);
#endif
#ifdef ATOMIC
					atom_add(&d_rhs_OSEM[local_ind], convert_ulong_sat(local_ele* *ax* TH));
#else
					atomicAdd_g_f(&d_rhs_OSEM[local_ind], (local_ele* *ax));
#endif
					local_ind += d_NN;
				}
			}
			else {
				local_ele *= *temp;
				for (uint kk = 0u; kk < d_N2; kk++) {
#ifdef ATOMIC
					atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
#else
					atomicAdd_g_f(&Summ[local_ind], local_ele);
#endif
					local_ind += d_NN;
				}
			}
		}
	}
	if (FP_bool) {
		*temp = 1. / *temp;
		if (d_attenuation_correction)
			* temp *= native_exp(jelppi);
		if (d_normalization == 1u)
			* temp *= d_norm;
		*temp *= global_factor;
	}
}