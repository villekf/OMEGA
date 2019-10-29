/**************************************************************************
* Functions for implementation 2 kernel files. No atomics. Used by both
* f32 and integer 64-bit atomic functions.
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
#define THR 0.001f

// Struct for boolean operators indicating whether a certain method is selected (OpenCL)
typedef struct _RecMethodsOpenCL {
	char MLEM, OSEM, MRAMLA, RAMLA, ROSEM, RBI, DRAMA, COSEM, ECOSEM, ACOSEM;
	char MRP, Quad, L, FMH, WeightedMean, TV, AD, APLS, TGV;
	char OSLMLEM, OSLOSEM, MBSREM, BSREM, ROSEMMAP, RBIMAP, OSLCOSEM;
} RecMethodsOpenCL;


// Denominator (forward projection)
inline void denominator(float local_ele, __local float* ax, const uint local_ind, const uint lid, const uint n_rekos,
	const uint d_N, const __global float* d_OSEM) {
	//uchar MAP = 0;
	uint yy = local_ind;
	for (uint kk = lid * n_rekos; kk < (lid + 1u) * n_rekos; kk++) {
		ax[kk] += (local_ele * d_OSEM[yy]);
		yy += d_N;
	}
}

// Nominator (backprojection) in MLEM
inline void nominator_mlem(__local float* ax, const float d_Sino, const float d_epps, const float temp, const uint randoms_correction, 
	const __global float* d_sc_ra, const uint idx, const uint lid, const uint n_rekos) {
	for (uint kk = lid * n_rekos; kk < (lid + 1u) * n_rekos; kk++) {
		if (ax[kk] <= 0.f)
			ax[kk] = d_epps;
		else
			ax[kk] *= temp;
		if (randoms_correction == 1u)
			ax[kk] += d_sc_ra[idx];
		ax[kk] = d_Sino / ax[kk];
	}
}

// Nominator (backprojection)
inline void nominator(__constant uchar* MethodList, __local float* ax, const float d_Sino, const float d_epsilon_mramla, const float d_epps, 
	const float temp, const uint randoms_correction, const __global float* d_sc_ra, const uint idx, const uint lid, const uint n_rekos) {
	uint oo = 0u;
	const float local_rand = d_sc_ra[idx];
	for (uint kk = lid * n_rekos; kk < (lid + 1u) * n_rekos; kk++) {
		if (ax[kk] <= 0.f)
			ax[kk] = d_epps;
		else
			ax[kk] *= temp;
		if (randoms_correction == 1u)
			ax[kk] += local_rand;
		if (MethodList[oo] != 1u)
			ax[kk] = d_Sino / ax[kk];
		else if (MethodList[oo] == 1u) { // MRAMLA/MBSREM
			if (randoms_correction == 1u) {
				if (ax[kk] <= d_epsilon_mramla && local_rand == 0.f && d_Sino > 0.f)
					ax[kk] = d_Sino / d_epsilon_mramla - (d_Sino / native_powr(d_epsilon_mramla, 2)) * (ax[kk] - d_epsilon_mramla);
				else
					ax[kk] = d_Sino / ax[kk];
			}
			else {
				if (ax[kk] <= d_epsilon_mramla && d_Sino > 0.f)
					ax[kk] = d_Sino / d_epsilon_mramla - (d_Sino / native_powr(d_epsilon_mramla, 2)) * (ax[kk] - d_epsilon_mramla);
				else
					ax[kk] = (d_Sino / ax[kk]);
			}
		}
		//ax[kk] = 0.f;
		oo++;
	}
}

 //Denominator (forward projection), orthogonal distance based ray tracer
inline void orth_distance_denominator(const int tempi, const uint Nx, const float y_diff, const float x_diff, const float y_center, __constant float* x_center, 
	const float kerroin, const float length_, float* temp, const uint temp_ijk, const float local_sino, __local float* ax, const uint d_Ny, const uint d_N, 
	const uint lid, const uint n_rekos, const uint im_dim, const __global float* d_OSEM) {

	const float diff = kerroin - x_diff * y_center;
	for (int uu = tempi; uu >= 0; uu--) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		*temp += local_ele;
		if (local_sino > 0.f) {
			denominator(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_OSEM);
		}
	}
	for (int uu = tempi + 1; uu < Nx; uu++) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		*temp += local_ele;
		if (local_sino > 0.f) {
			denominator(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_OSEM);
		}
	}
}


// Denominator (forward projection), orthogonal distance based ray tracer
inline void orth_distance_denominator_3D(const int tempi, const uint Nx, const uint Nz, const float y_diff, const float x_diff, const float z_diff,
	const float y_center, __constant float* x_center, __constant float* z_center, const float crystal_size_z, float* temp, const uint temp_ijk, 
	const float local_sino, __local float* ax, const uint d_Ny, const uint d_N, const int tempk, const uint Nxy, const float xs, const float ys, const float zs,
	const int dec, const uint lid, const uint n_rekos, const uint im_dim,
	const __global float* d_OSEM) {

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
			*temp += local_ele;
			if (local_sino > 0.f) {
				denominator(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_OSEM);
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
				*temp += local_ele;
				if (local_sino > 0.f) {
					denominator(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_OSEM);
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
			*temp += local_ele;
			if (local_sino > 0.f) {
				denominator(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_OSEM);
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
				*temp += local_ele;
				if (local_sino > 0.f) {
					denominator(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_OSEM);
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
inline void orth_distance_denominator_perpendicular(const float diff1, __constant float* center1, const float kerroin, const float length_, float* temp, 
	const uint d_attenuation_correction, const uint d_normalization, __local float* ax, const float d_b, const float d, const float d_d1, const uint d_N1, 
	const uint d_N2, const uint z_loop, const __global float* d_atten, const __global float* d_norm,  const uint idx, const float local_sino, 
	const uint d_N, const uint d_NN, const uint lid, const uint n_rekos, const uint im_dim, const __global float* d_OSEM) {

	const uint zz = z_loop * d_N2 * d_N1;
	const int apu = perpendicular_start(d_b, d, d_d1, d_N1);
	float jelppi = 0.f;
	for (int uu = apu; uu >= 0; uu--) {
		float local_ele = compute_element_orth(kerroin, diff1, center1[uu], length_);
		if (local_ele <= THR)
			break;
		*temp += (local_ele * d_N2);
		uint local_ind = convert_uint_sat(uu) * d_N + zz;
		for (uint kk = 0u; kk < d_N2; kk++) {
			if (d_attenuation_correction && uu == apu)
				jelppi += (d_d1 * -d_atten[local_ind]);
			if (local_sino > 0.f) {
				denominator(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_OSEM);
			}
			local_ind += d_NN;
		}
	}
	for (uint uu = convert_uint_sat(apu) + 1u; uu < d_N1; uu++) {
		float local_ele = compute_element_orth(kerroin, diff1, center1[uu], length_);
		if (local_ele <= THR)
			break;
		*temp += (local_ele * d_N2);
		uint local_ind = uu * d_N + zz;
		for (uint kk = 0u; kk < d_N2; kk++) {
			if (local_sino > 0.f) {
				denominator(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_OSEM);
			}
			local_ind += d_NN;
		}
	}
	*temp = 1.f / *temp;
	if (d_attenuation_correction == 1u)
		* temp *= jelppi;
	if (d_normalization == 1u)
		* temp *= d_norm[idx];
}


inline void orth_distance_denominator_perpendicular_3D(__constant float* center1, const float center2, __constant float* z_center, float* temp, 
	const uint d_attenuation_correction, const uint d_normalization, __local float* ax, const float d_b, const float d, const float d_d1, const uint d_N1, 
	const uint d_N2, const uint z_loop,	const __global float* d_atten, const __global float* d_norm, const uint idx, const float local_sino, const uint d_N,
	const uint d_NN, const float xs, const float ys, const float zs, const float xl, const float yl, const float zl, const float crystal_size_z, const uint Nyx,
	const uint Nz, const uint lid, const uint n_rekos, const uint im_dim, const __global float* d_OSEM) {

	const uint apu = perpendicular_start(d_b, d, d_d1, d_N1);
	float jelppi = 0.f;
	for (int zz = convert_int_sat(z_loop); zz >= 0; zz--) {
		for (int uu = convert_int_sat(apu); uu >= 0; uu--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			*temp += (local_ele * d_N2);
			for (uint kk = 0u; kk < d_N2; kk++) {
				if (d_attenuation_correction && zz == convert_int_sat(z_loop) && uu == convert_int_sat(apu))
					jelppi += (d_d1 * -d_atten[local_ind]);
				if (local_sino > 0.f) {
					denominator(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_OSEM);
				}
				local_ind += d_NN;
			}
		}
		for (uint uu = apu + 1; uu < d_N1; uu++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			*temp += (local_ele * d_N2);
			for (uint kk = 0u; kk < d_N2; kk++) {
				if (d_attenuation_correction && zz == convert_int_sat(z_loop) && uu == convert_int_sat(apu))
					jelppi += (d_d1 * -d_atten[local_ind]);
				if (local_sino > 0.f) {
					denominator(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_OSEM);
				}
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
			*temp += (local_ele * d_N2);
			for (uint kk = 0u; kk < d_N2; kk++) {
				if (d_attenuation_correction && zz == convert_int_sat(z_loop) && uu == convert_int_sat(apu))
					jelppi += (d_d1 * -d_atten[local_ind]);
				if (local_sino > 0.f) {
					denominator(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_OSEM);
				}
				local_ind += d_NN;
			}
		}
		for (uint uu = apu + 1; uu < d_N1; uu++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			*temp += (local_ele * d_N2);
			for (uint kk = 0u; kk < d_N2; kk++) {
				if (d_attenuation_correction && zz == convert_int_sat(z_loop) && uu == convert_int_sat(apu))
					jelppi += (d_d1 * -d_atten[local_ind]);
				if (local_sino > 0.f) {
					denominator(local_ele, ax, local_ind, lid, n_rekos, im_dim, d_OSEM);
				}
				local_ind += d_NN;
			}
		}
	}
	*temp = 1.f / *temp;
	if (d_attenuation_correction == 1u)
		* temp *= jelppi;
	if (d_normalization == 1u)
		* temp *= d_norm[idx];
}

// Nominator (backprojection), COSEM
inline void nominator_cosem(float* axCOSEM, const float local_sino, const float d_epps, const float temp, const uint d_randoms, const __global float* d_sc_ra,
	const uint idx) {
	if (*axCOSEM <= 0.f)
		*axCOSEM = d_epps;
	else
		*axCOSEM *= temp;
	if (d_randoms == 1u)
		*axCOSEM += d_sc_ra[idx];
	*axCOSEM = local_sino / *axCOSEM;
}