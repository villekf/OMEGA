#pragma once
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
#include "general_opencl_functions.h"

#define THR 0.001f
#define TH 100000000000.f

// Denominator (forward projection), multi-GPU version, orthogonal distance based ray tracer
inline void orth_distance_multi_a64(const int tempi, const uint Nx, const float y_diff, const float x_diff, const float y_center,
	__constant float* x_center, const float kerroin, const float length_, float* temp, const uint temp_ijk,
	const float local_sino, float* ax, const uint d_Ny, const uint d_N, const uchar no_norm, const bool RHS, const bool SUMMA, const __global float* d_OSEM, __global ulong* Summ, __global ulong* d_rhs_OSEM) {

	const float diff = kerroin - x_diff * y_center;
	for (int uu = tempi; uu >= 0; uu--) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		if (RHS) {
			local_ele *= *temp;
			atom_add(&d_rhs_OSEM[local_ind], convert_ulong_sat(local_ele * *ax * TH));
			if (no_norm == 0u)
				atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
		}
		else if (SUMMA) {
			local_ele *= *temp;
			atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
		}
		else {
			*temp += local_ele;
			if (local_sino > 0.f) {
				denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
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
			atom_add(&d_rhs_OSEM[local_ind], convert_ulong_sat(local_ele * *ax * TH));
			if (no_norm == 0u)
				atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
		}
		else if (SUMMA) {
			local_ele *= *temp;
			atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
		}
		else {
			*temp += local_ele;
			if (local_sino > 0.f) {
				denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
			}
		}
	}
}


// Denominator (forward projection), multi-GPU version, orthogonal distance based ray tracer
inline void orth_distance_multi_3D_a64(const int tempi, const uint Nx, const uint Nz, const float y_diff, const float x_diff, const float z_diff,
	const float y_center, __constant float* x_center, __constant float* z_center, const float crystal_size_z, float* temp, const uint temp_ijk,
	const float local_sino, float* ax, const uint d_Ny, const uint d_N, const int tempk, const uint Nxy, const float xs, const float ys, const float zs,
	const int dec, const uchar no_norm, const bool RHS, const bool SUMMA, const __global float* d_OSEM, __global ulong* Summ, __global ulong* d_rhs_OSEM) {

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
				atom_add(&d_rhs_OSEM[local_ind], convert_ulong_sat(local_ele * *ax * TH));
				if (no_norm == 0u)
					atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
			}
			else if (SUMMA) {
				local_ele *= *temp;
				atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
			}
			else {
				*temp += local_ele;
				if (local_sino > 0.f) {
					denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
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
					atom_add(&d_rhs_OSEM[local_ind], convert_ulong_sat(local_ele * *ax * TH));
					if (no_norm == 0u)
						atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
				}
				else if (SUMMA) {
					local_ele *= *temp;
					atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
				}
				else {
					*temp += local_ele;
					if (local_sino > 0.f) {
						denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
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
				atom_add(&d_rhs_OSEM[local_ind], convert_ulong_sat(local_ele * *ax * TH));
				if (no_norm == 0u)
					atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
			}
			else if (SUMMA) {
				local_ele *= *temp;
				atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
			}
			else {
				*temp += local_ele;
				if (local_sino > 0.f) {
					denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
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
					atom_add(&d_rhs_OSEM[local_ind], convert_ulong_sat(local_ele * *ax * TH));
					if (no_norm == 0u)
						atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
				}
				else if (SUMMA) {
					local_ele *= *temp;
					atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
				}
				else {
					*temp += local_ele;
					if (local_sino > 0.f) {
						denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
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


//// Denominator, forward/back projection, orthogonal distance based ray tracer
//inline void orth_distance_bpfp_a64(const int tempi, const uint Nx, const float y_diff, const float x_diff, const float y_center,
//	__constant float* x_center, const float kerroin, const float length_, float* temp, const uint temp_ijk,
//	const uint d_Ny, const uint d_N, const int tempj, const uchar fp, float* d_output_local, const __global float* d_f, __global float* d_Summ,
//	__global ulong* bp_output, const uchar no_norm, float d_rhs_local, const bool FP_bool) {
//
//	const float diff = kerroin - x_diff * y_center;
//	for (int uu = tempi; uu >= 0; uu--) {
//		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
//		if (local_ele <= THR)
//			break;
//		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
//		if (FP_bool) {
//			*temp += local_ele;
//			if (fp) {
//				*d_output_local += (d_f[local_ind] * local_ele);
//			}
//		}
//		else {
//			local_ele *= *temp;
//			atom_add(&bp_output[local_ind], convert_ulong_sat(local_ele * d_rhs_local * TH));
//			if (!no_norm)
//				atom_add(&d_Summ[local_ind], convert_ulong_sat(local_ele * TH));
//		}
//	}
//	for (int uu = tempi + 1; uu < Nx; uu++) {
//		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
//		if (local_ele <= THR)
//			break;
//		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
//		if (FP_bool) {
//			*temp += local_ele;
//			if (fp) {
//				*d_output_local += (d_f[local_ind] * local_ele);
//			}
//		}
//		else {
//			local_ele *= *temp;
//			atom_add(&bp_output[local_ind], convert_ulong_sat(local_ele * d_rhs_local * TH));
//			if (!no_norm)
//				atom_add(&d_Summ[local_ind], convert_ulong_sat(local_ele * TH));
//		}
//	}
//}


//// Denominator, forward/back projection, orthogonal distance based ray tracer
//inline void orth_distance_bpfp_3D_a64(const int tempi, const uint Nx, const uint Nz, const float y_diff, const float x_diff, const float z_diff,
//	const float y_center, __constant float* x_center, __constant float* z_center, const float crystal_size_z, float* temp, const uint temp_ijk,
//	const uint d_Ny, const uint d_N, const int tempk, const uint Nxy, const float xs, const float ys, const float zs, const int dec, const uchar fp, float* d_output_local, const __global float* d_f, __global float* d_Summ,
//	__global ulong* bp_output, const uchar no_norm, float d_rhs_local, const bool FP_bool) {
//
//	bool loppu1 = true;
//	bool loppu2 = true;
//	uchar loppu3 = 0u;
//	uchar loppu4 = 0u;
//	bool pass1 = true;
//	bool pass = false;
//	int alku = tempk;
//	int alku2 = tempk;
//	int alku3 = tempk;
//	int alku4 = tempk + 1;
//	for (int uu = tempi; uu >= 0; uu--) {
//		//for (int zz = tempk; zz >= 0; zz--) {
//		for (int zz = alku; zz >= 0; zz--) {
//			float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
//			if (local_ele <= THR) {
//				if (loppu1 && zz >= alku3 - dec && !(pass && alku3 == tempk)) {
//					continue;
//				}
//				else {
//					if (loppu1 && loppu3 == 1u && loppu4 == 0u)
//						pass1 = false;
//					break;
//				}
//			}
//			const uint local_ind = compute_ind_orth_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
//			if (FP_bool) {
//				*temp += local_ele;
//				if (fp) {
//					*d_output_local += (d_f[local_ind] * local_ele);
//				}
//			}
//			else {
//				local_ele *= *temp;
//				atom_add(&bp_output[local_ind], convert_ulong_sat(local_ele * d_rhs_local * TH));
//				if (!no_norm)
//					atom_add(&d_Summ[local_ind], convert_ulong_sat(local_ele * TH));
//			}
//			if (loppu1) {
//				alku = zz;
//				loppu1 = false;
//				loppu3 = 1u;
//				if (zz == tempk || loppu4 == 1u)
//					pass = true;
//			}
//			alku3 = zz;
//		}
//		if (pass) {
//			for (int zz = alku2 + 1; zz < Nz; zz++) {
//				//for (uint zz = convert_uint_sat(tempk) + 1u; zz < Nz; zz++) {
//				float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
//				if (local_ele <= THR) {
//					if (loppu2 && zz <= alku4 + dec) {
//						continue;
//					}
//					else {
//						if (loppu2 && loppu4 == 1u) {
//							loppu4 = 0u;
//							pass = false;
//						}
//						break;
//					}
//				}
//				const uint local_ind = compute_ind_orth_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
//				if (FP_bool) {
//					*temp += local_ele;
//					if (fp) {
//						*d_output_local += (d_f[local_ind] * local_ele);
//					}
//				}
//				else {
//					local_ele *= *temp;
//					atom_add(&bp_output[local_ind], convert_ulong_sat(local_ele * d_rhs_local * TH));
//					if (!no_norm)
//						atom_add(&d_Summ[local_ind], convert_ulong_sat(local_ele * TH));
//				}
//				if (loppu2) {
//					alku2 = zz - 1;
//					loppu2 = false;
//					loppu4 = 1u;
//					pass1 = true;
//				}
//				alku4 = zz;
//			}
//		}
//		loppu1 = true;
//		loppu2 = true;
//		if (!pass1) {
//			break;
//		}
//	}
//	loppu3 = 0u;
//	pass1 = true;
//	pass = false;
//	loppu4 = 0u;
//	alku2 = tempk;
//	alku = tempk;
//	alku3 = tempk;
//	alku4 = tempk + 1;
//	for (uint uu = convert_uint_sat(tempi) + 1u; uu < Nx; uu++) {
//		for (int zz = alku; zz >= 0; zz--) {
//			//for (int zz = tempk; zz >= 0; zz--) {
//			float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
//			if (local_ele <= THR) {
//				if (loppu1 && zz >= alku3 - dec && !(pass && alku3 == tempk)) {
//					continue;
//				}
//				else {
//					if (loppu1 && loppu3 == 1u && loppu4 == 0u)
//						pass1 = false;
//					break;
//				}
//			}
//			const uint local_ind = compute_ind_orth_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
//			if (FP_bool) {
//				*temp += local_ele;
//				if (fp) {
//					*d_output_local += (d_f[local_ind] * local_ele);
//				}
//			}
//			else {
//				local_ele *= *temp;
//				atom_add(&bp_output[local_ind], convert_ulong_sat(local_ele * d_rhs_local * TH));
//				if (!no_norm)
//					atom_add(&d_Summ[local_ind], convert_ulong_sat(local_ele * TH));
//			}
//			if (loppu1) {
//				alku = zz;
//				loppu1 = false;
//				loppu3 = 1u;
//				if (zz == tempk || loppu4 == 1u)
//					pass = true;
//			}
//			alku3 = zz;
//		}
//		if (pass) {
//			for (int zz = alku2 + 1; zz < Nz; zz++) {
//				//for (uint zz = convert_uint_sat(tempk) + 1u; zz < Nz; zz++) {
//				float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
//				if (local_ele <= THR) {
//					if (loppu2 && zz <= alku4 + dec) {
//						continue;
//					}
//					else {
//						if (loppu2 && loppu4 == 1u) {
//							loppu4 = 0u;
//							pass = false;
//						}
//						break;
//					}
//				}
//				const uint local_ind = compute_ind_orth_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
//				if (FP_bool) {
//					*temp += local_ele;
//					if (fp) {
//						*d_output_local += (d_f[local_ind] * local_ele);
//					}
//				}
//				else {
//					local_ele *= *temp;
//					atom_add(&bp_output[local_ind], convert_ulong_sat(local_ele * d_rhs_local * TH));
//					if (!no_norm)
//						atom_add(&d_Summ[local_ind], convert_ulong_sat(local_ele * TH));
//				}
//				if (loppu2) {
//					alku2 = zz - 1;
//					loppu2 = false;
//					loppu4 = 1u;
//					pass1 = true;
//				}
//				alku4 = zz;
//			}
//		}
//		loppu1 = true;
//		loppu2 = true;
//		if (!pass1) {
//			break;
//		}
//	}
//}

// Denominator (forward projection), orthogonal distance based ray tracer, multi-GPU
inline void orth_distance_perpendicular_multi_a64(const float diff1, __constant float* center1, const float kerroin,
	const float length_, float* temp, const uint d_attenuation_correction, const uint d_normalization, float* ax, const float d_b, const float d, const float d_d1,
	const uint d_N1, const uint d_N2, const uint z_loop, const __global float* d_atten, const float d_norm, const float local_sino, const uint d_N, const uint d_NN,
	const __global float* d_OSEM, const uchar no_norm, __global ulong* d_rhs_OSEM, __global ulong* Summ, const bool FP, const bool RHS) {

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
					denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
				}
				local_ind += d_NN;
			}
		}
		else if (RHS) {
			local_ele *= *temp;
			for (uint kk = 0u; kk < d_N2; kk++) {
				if (no_norm == 0u)
					atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
				atom_add(&d_rhs_OSEM[local_ind], convert_ulong_sat(local_ele * *ax * TH));
				local_ind += d_NN;
			}
		}
		else {
			local_ele *= *temp;
			for (uint kk = 0u; kk < d_N2; kk++) {
				atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
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
					denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
				}
				local_ind += d_NN;
			}
		}
		else if (RHS) {
			local_ele *= *temp;
			for (uint kk = 0u; kk < d_N2; kk++) {
				if (no_norm == 0u)
					atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
				atom_add(&d_rhs_OSEM[local_ind], convert_ulong_sat(local_ele * *ax * TH));
				local_ind += d_NN;
			}
		}
		else {
			local_ele *= *temp;
			for (uint kk = 0u; kk < d_N2; kk++) {
				atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
				local_ind += d_NN;
			}
		}
	}
	if (FP) {
		*temp = 1.f / *temp;
		if (d_attenuation_correction)
			* temp *= jelppi;
		if (d_normalization == 1u)
			* temp *= d_norm;
	}
}

//// Denominator (forward projection), orthogonal distance based ray tracer, multi-GPU
//inline void orth_distance_perpendicular_bpfp_a64(const float diff1, __constant float* center1, const float kerroin,
//	const float length_, float* temp, const uint d_attenuation_correction, const uint d_normalization, const float d_b, const float d, const float d_d1,
//	const uint d_N1, const uint d_N2, const uint z_loop, const __global float* d_atten, const float d_norm, const uint d_N, const uint d_NN,
//	const __global float* d_rhs, const uchar fp, float* d_output_local, __global ulong* bp_output, const __global float* d_f, __global ulong* Summ, const uchar no_norm, const float d_rhs_local, const bool FP_bool) {
//
//	const uint zz = z_loop * d_N2 * d_N1;
//	const int apu = perpendicular_start(d_b, d, d_d1, d_N1);
//	float jelppi = 0.f;
//	for (int uu = apu; uu >= 0; uu--) {
//		float local_ele = compute_element_orth(kerroin, diff1, center1[uu], length_);
//		if (local_ele <= THR)
//			break;
//		uint local_ind = convert_uint_sat(uu) * d_N + zz;
//		if (FP_bool) {
//			*temp += (local_ele * d_N2);
//			if ((d_attenuation_correction && uu == apu) || fp) {
//				for (uint kk = 0u; kk < d_N2; kk++) {
//					if (d_attenuation_correction && uu == apu)
//						jelppi += (d_d1 * -d_atten[local_ind]);
//					if (fp)
//						* d_output_local += (d_f[local_ind] * local_ele);
//					local_ind += d_NN;
//				}
//			}
//		}
//		else {
//			local_ele *= *temp;
//			for (uint kk = 0u; kk < d_N2; kk++) {
//				atom_add(&bp_output[local_ind], convert_ulong_sat(local_ele * d_rhs_local * TH));
//				if (!no_norm)
//					atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
//				local_ind += d_NN;
//			}
//		}
//	}
//	for (uint uu = convert_uint_sat(apu) + 1u; uu < d_N1; uu++) {
//		float local_ele = compute_element_orth(kerroin, diff1, center1[uu], length_);
//		if (local_ele <= THR)
//			break;
//		uint local_ind = uu * d_N + zz;
//		if (FP_bool) {
//			*temp += (local_ele * d_N2);
//			if (fp) {
//				for (uint kk = 0u; kk < d_N2; kk++) {
//					*d_output_local += (d_f[local_ind] * local_ele);
//					local_ind += d_NN;
//				}
//			}
//		}
//		else {
//			local_ele *= *temp;
//			for (uint kk = 0u; kk < d_N2; kk++) {
//				atom_add(&bp_output[local_ind], convert_ulong_sat(local_ele * d_rhs_local * TH));
//				if (!no_norm)
//					atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
//				local_ind += d_NN;
//			}
//		}
//	}
//	if (FP_bool) {
//		*temp = 1.f / *temp;
//		if (d_attenuation_correction)
//			* temp *= jelppi;
//		if (d_normalization == 1u)
//			* temp *= d_norm;
//	}
//}

//// Denominator (forward projection), orthogonal distance based ray tracer, multi-GPU
//inline void orth_distance_perpendicular_bpfp_3D_a64(__constant float* center1, const float center2, __constant float* z_center,
//	float* temp, const uint d_attenuation_correction, const uint d_normalization, const float d_b, const float d, const float d_d1,
//	const uint d_N1, const uint d_N2, const uint z_loop, const __global float* d_atten, const float d_norm, const uint d_N, const uint d_NN, const __global float* d_rhs, const uchar fp, float* d_output_local,
//	const float xs, const float ys, const float zs, const float xl, const float yl, const float zl, const float crystal_size_z, const uint Nyx, const uint Nz,
//	__global ulong* bp_output, const __global float* d_f, __global ulong* Summ, const uchar no_norm, const float d_rhs_local, const bool FP_bool) {
//
//	const uint apu = perpendicular_start(d_b, d, d_d1, d_N1);
//	float jelppi = 0.f;
//	for (int zz = convert_int_sat(z_loop); zz >= 0; zz--) {
//		for (int uu = convert_int_sat(apu); uu >= 0; uu--) {
//			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
//			if (local_ele <= THR)
//				break;
//			uint local_ind = uu * d_N + zz * Nyx;
//			if (FP_bool) {
//				*temp += (local_ele * d_N2);
//				if ((d_attenuation_correction && zz == convert_int_sat(z_loop) && uu == convert_int_sat(apu)) || fp) {
//					for (uint kk = 0u; kk < d_N2; kk++) {
//						if (d_attenuation_correction && zz == convert_int_sat(z_loop) && uu == convert_int_sat(apu))
//							jelppi += (d_d1 * -d_atten[local_ind]);
//						if (fp)
//							* d_output_local += (d_f[local_ind] * local_ele);
//						local_ind += d_NN;
//					}
//				}
//			}
//			else {
//				local_ele *= *temp;
//				for (uint kk = 0u; kk < d_N2; kk++) {
//					atom_add(&bp_output[local_ind], convert_ulong_sat(local_ele * d_rhs_local * TH));
//					if (!no_norm)
//						atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
//					local_ind += d_NN;
//				}
//			}
//		}
//		for (uint uu = convert_uint_sat(apu) + 1u; uu < d_N1; uu++) {
//			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
//			if (local_ele <= THR)
//				break;
//			uint local_ind = uu * d_N + zz * Nyx;
//			if (FP_bool) {
//				*temp += (local_ele * d_N2);
//				if (fp) {
//					for (uint kk = 0u; kk < d_N2; kk++) {
//						*d_output_local += (d_f[local_ind] * local_ele);
//						local_ind += d_NN;
//					}
//				}
//			}
//			else {
//				local_ele *= *temp;
//				for (uint kk = 0u; kk < d_N2; kk++) {
//					atom_add(&bp_output[local_ind], convert_ulong_sat(local_ele * d_rhs_local * TH));
//					if (!no_norm)
//						atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
//					local_ind += d_NN;
//				}
//			}
//		}
//	}
//	for (uint zz = z_loop + 1u; zz < Nz; zz++) {
//		for (int uu = convert_int_sat(apu); uu >= 0; uu--) {
//			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
//			if (local_ele <= THR)
//				break;
//			uint local_ind = uu * d_N + zz * Nyx;
//			if (FP_bool) {
//				*temp += (local_ele * d_N2);
//				if (fp) {
//					for (uint kk = 0u; kk < d_N2; kk++) {
//						*d_output_local += (d_f[local_ind] * local_ele);
//						local_ind += d_NN;
//					}
//				}
//			}
//			else {
//				local_ele *= *temp;
//				for (uint kk = 0u; kk < d_N2; kk++) {
//					atom_add(&bp_output[local_ind], convert_ulong_sat(local_ele * d_rhs_local * TH));
//					if (!no_norm)
//						atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
//					local_ind += d_NN;
//				}
//			}
//		}
//		for (uint uu = convert_uint_sat(apu) + 1u; uu < d_N1; uu++) {
//			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
//			if (local_ele <= THR)
//				break;
//			uint local_ind = uu * d_N + zz * Nyx;
//			if (FP_bool) {
//				*temp += (local_ele * d_N2);
//				if (fp) {
//					for (uint kk = 0u; kk < d_N2; kk++) {
//						*d_output_local += (d_f[local_ind] * local_ele);
//						local_ind += d_NN;
//					}
//				}
//			}
//			else {
//				local_ele *= *temp;
//				for (uint kk = 0u; kk < d_N2; kk++) {
//					atom_add(&bp_output[local_ind], convert_ulong_sat(local_ele * d_rhs_local * TH));
//					if (!no_norm)
//						atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
//					local_ind += d_NN;
//				}
//			}
//		}
//	}
//	if (FP_bool) {
//		*temp = 1.f / *temp;
//		if (d_attenuation_correction)
//			* temp *= jelppi;
//		if (d_normalization == 1u)
//			* temp *= d_norm;
//	}
//}

inline void orth_distance_perpendicular_multi_3D_a64(__constant float* center1, const float center2, __constant float* z_center,
	float* temp, const uint d_attenuation_correction, const uint d_normalization, float* ax, const float d_b, const float d, const float d_d1,
	const uint d_N1, const uint d_N2, const uint z_loop, const __global float* d_atten, const float d_norm, const float local_sino, const uint d_N, const uint d_NN, const __global float* d_OSEM,
	const float xs, const float ys, const float zs, const float xl, const float yl, const float zl, const float crystal_size_z, const uint Nyx, const uint Nz,
	const uchar no_norm, __global ulong* Summ, __global ulong* d_rhs_OSEM, const bool FP, const bool RHS) {

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
						denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
					}
					local_ind += d_NN;
				}
			}
			else if (RHS) {
				local_ele *= *temp;
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (no_norm == 0u)
						atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
					atom_add(&d_rhs_OSEM[local_ind], convert_ulong_sat(local_ele * *ax * TH));
					local_ind += d_NN;
				}
			}
			else {
				local_ele *= *temp;
				for (uint kk = 0u; kk < d_N2; kk++) {
					atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
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
						denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
					}
					local_ind += d_NN;
				}
			}
			else if (RHS) {
				local_ele *= *temp;
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (no_norm == 0u)
						atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
					atom_add(&d_rhs_OSEM[local_ind], convert_ulong_sat(local_ele * *ax * TH));
					local_ind += d_NN;
				}
			}
			else {
				local_ele *= *temp;
				for (uint kk = 0u; kk < d_N2; kk++) {
					atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
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
						denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
					}
					local_ind += d_NN;
				}
			}
			else if (RHS) {
				local_ele *= *temp;
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (no_norm == 0u)
						atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
					atom_add(&d_rhs_OSEM[local_ind], convert_ulong_sat(local_ele * *ax * TH));
					local_ind += d_NN;
				}
			}
			else {
				local_ele *= *temp;
				for (uint kk = 0u; kk < d_N2; kk++) {
					atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
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
						denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
					}
					local_ind += d_NN;
				}
			}
			else if (RHS) {
				local_ele *= *temp;
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (no_norm == 0u)
						atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
					atom_add(&d_rhs_OSEM[local_ind], convert_ulong_sat(local_ele * *ax * TH));
					local_ind += d_NN;
				}
			}
			else {
				local_ele *= *temp;
				for (uint kk = 0u; kk < d_N2; kk++) {
					atom_add(&Summ[local_ind], convert_ulong_sat(local_ele * TH));
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
			* temp *= d_norm;
	}
}