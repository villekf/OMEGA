/**************************************************************************
* Special functions for the volume-based ray tracer.
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
#pragma once

__device__ void vol_perpendicular_multi_3D(const float* center1, const float center2, const float* z_center,
	float* temp, const unsigned int d_attenuation_correction, const unsigned int d_normalization, float* ax, const float d_b, const float d, const float d_d1,
	const unsigned int d_N1, const unsigned int d_N2, const unsigned int z_loop, const float* d_atten, const float d_norm, const float local_sino, const unsigned int d_N,
	const unsigned int d_NN, const float* d_OSEM, const float xs, const float ys, const float zs, const float xl, const float yl, const float zl,
	const float crystal_size_z, const unsigned int Nyx, const unsigned int Nz, const unsigned char no_norm, CAST* Summ, CAST* d_rhs_OSEM,
	const bool FP, const bool RHS, const float global_factor, const float bmin, const float bmax, const float Vmax, const float* V) {

	//const unsigned int zz = z_loop * d_N2 * d_N1;
	const unsigned int apu = perpendicular_start(d_b, d, d_d1, d_N1);
	float jelppi = 0.f;
	for (int zz = (int)(z_loop); zz >= 0; zz--) {
		for (int uu = (int)(apu); uu >= 0; uu--) {
			float local_ele = compute_element_vol_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele >= bmax)
				break;
			unsigned int local_ind = uu * d_N + zz * Nyx;
			if (local_ele < bmin)
				local_ele = Vmax;
			else
				local_ele = V[__float2uint_rn((local_ele - bmin) * CC)];
			if (FP) {
				*temp += (local_ele * d_N2);
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (d_attenuation_correction && zz == (int)(z_loop) && uu == (int)(apu))
						jelppi += (d_d1 * -d_atten[local_ind]);
					if (local_sino > 0.f) {
						denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
					}
					local_ind += d_NN;
				}
			}
			else if (RHS) {
				local_ele *= *temp;
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (no_norm == 0u)
#ifdef ATOMIC
						atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
						atomicAdd(&Summ[local_ind], local_ele);
#endif
#ifdef ATOMIC
					atomicAdd(&d_rhs_OSEM[local_ind], __float2ull_rn(local_ele * *ax * TH));
#else
					atomicAdd(&d_rhs_OSEM[local_ind], (local_ele * *ax));
#endif
					local_ind += d_NN;
				}
			}
			else {
				local_ele *= *temp;
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					atomicAdd(&Summ[local_ind], local_ele);
					local_ind += d_NN;
				}
			}
		}
		for (unsigned int uu = apu + 1; uu < d_N1; uu++) {
			float local_ele = compute_element_vol_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele >= bmax)
				break;
			if (local_ele < bmin)
				local_ele = Vmax;
			else
				local_ele = V[__float2uint_rn((local_ele - bmin) * CC)];
			unsigned int local_ind = uu * d_N + zz * Nyx;
			if (FP) {
				*temp += (local_ele * d_N2);
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (local_sino > 0.f) {
						denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
					}
					local_ind += d_NN;
				}
			}
			else if (RHS) {
				local_ele *= *temp;
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (no_norm == 0u)
#ifdef ATOMIC
						atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
						atomicAdd(&Summ[local_ind], local_ele);
#endif
#ifdef ATOMIC
					atomicAdd(&d_rhs_OSEM[local_ind], __float2ull_rn(local_ele * *ax * TH));
#else
					atomicAdd(&d_rhs_OSEM[local_ind], (local_ele * *ax));
#endif
					local_ind += d_NN;
				}
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
	for (unsigned int zz = z_loop + 1u; zz < Nz; zz++) {
		for (int uu = (int)(apu); uu >= 0; uu--) {
			float local_ele = compute_element_vol_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele >= bmax)
				break;
			if (local_ele < bmin)
				local_ele = Vmax;
			else
				local_ele = V[__float2uint_rn((local_ele - bmin) * CC)];
			unsigned int local_ind = uu * d_N + zz * Nyx;
			if (FP) {
				*temp += (local_ele * d_N2);
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (local_sino > 0.f) {
						denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
					}
					local_ind += d_NN;
				}
			}
			else if (RHS) {
				local_ele *= *temp;
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (no_norm == 0u)
#ifdef ATOMIC
						atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
						atomicAdd(&Summ[local_ind], local_ele);
#endif
#ifdef ATOMIC
					atomicAdd(&d_rhs_OSEM[local_ind], __float2ull_rn(local_ele * *ax * TH));
#else
					atomicAdd(&d_rhs_OSEM[local_ind], (local_ele * *ax));
#endif
					local_ind += d_NN;
				}
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
			float local_ele = compute_element_vol_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele >= bmax)
				break;
			if (local_ele < bmin)
				local_ele = Vmax;
			else
				local_ele = V[__float2uint_rn((local_ele - bmin) * CC)];
			unsigned int local_ind = uu * d_N + zz * Nyx;
			if (FP) {
				*temp += (local_ele * d_N2);
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (local_sino > 0.f) {
						denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
					}
					local_ind += d_NN;
				}
			}
			else if (RHS) {
				local_ele *= *temp;
				for (unsigned int kk = 0u; kk < d_N2; kk++) {
					if (no_norm == 0u)
#ifdef ATOMIC
						atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
						atomicAdd(&Summ[local_ind], local_ele);
#endif
#ifdef ATOMIC
					atomicAdd(&d_rhs_OSEM[local_ind], __float2ull_rn(local_ele * *ax * TH));
#else
					atomicAdd(&d_rhs_OSEM[local_ind], (local_ele * *ax));
#endif
					local_ind += d_NN;
				}
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
	if (FP) {
		*temp = 1. / *temp;
		if (d_attenuation_correction)
			*temp *= native_exp(jelppi);
		if (d_normalization == 1u)
			*temp *= d_norm;
		*temp *= global_factor;
	}
}