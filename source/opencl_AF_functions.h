#pragma once
#define THR 0.001f

// Struct for boolean operators indicating whether a certain method is selected (OpenCL)
typedef struct _RecMethodsOpenCL {
	char MLEM, OSEM, MRAMLA, RAMLA, ROSEM, RBI, DRAMA, COSEM, ECOSEM, ACOSEM;
	char MRP, Quad, L, FMH, WeightedMean, TV, AD, APLS, TGV;
	char OSLMLEM, OSLOSEM, MBSREM, BSREM, ROSEMMAP, RBIMAP, OSLCOSEM;
} RecMethodsOpenCL;

// This function was taken from: https://streamhpc.com/blog/2016-02-09/atomic-operations-for-floats-in-opencl-improved/
// Computes the atomic_add for floats
inline void atomicAdd_g_f(volatile __global float *addr, float val) {
	union {
		unsigned int u32;
		float        f32;
	} next, expected, current;
	current.f32 = *addr;
	do {
		expected.f32 = current.f32;
		next.f32 = expected.f32 + val;
		current.u32 = atomic_cmpxchg((volatile __global unsigned int *)addr, expected.u32, next.u32);
	} while (current.u32 != expected.u32);
}

// Compute the Euclidean norm of a vector
inline float e_norm(const float x, const float y, const float z) {
	return native_sqrt(x * x + y * y + z * z);
}

inline float compute_element_orth_3D(const float xs, const float ys, const float zs, const float xl, const float yl, const float zl, const float crystal_size_z,
	const float xp, const float yp, const float zp) {

	float x1, y1, z1, x0, y0, z0;

	x0 = xp - xs;
	y0 = yp - ys;
	z0 = zp - zs;

	// Cross product
	x1 = yl * z0 - zl * y0;
	y1 = zl * x0 - xl * z0;
	z1 = xl * y0 - yl * x0;

	const float normi = e_norm(x1, y1, z1);

	return (1.f - normi / crystal_size_z);
}

inline uint compute_ind_orth_mfree_3D(const uint tempi, const uint tempijk, const int tempk, const uint d_N, const uint Nyx) {
	uint local_ind = tempi * d_N + tempijk + convert_uint_sat(tempk) * Nyx;
	return local_ind;
}

// Denominator (forward projection) in MLEM
inline void denominator_mlem(float local_ele, const RecMethodsOpenCL MethodList, float ax[], const uint local_ind, const __global float* d_MLEM,
	const __global float* d_MRP_MLEM, const __global float* d_Quad_MLEM, const __global float* d_L_MLEM, const __global float* d_FMH_MLEM, const __global float* d_WeightedMean_MLEM,
	const __global float* d_TV_MLEM, const __global float* d_AD_MLEM, const __global float* d_APLS_MLEM, const __global float* d_TGV_MLEM) {
	if (MethodList.MLEM)
		ax[0] += (local_ele * d_MLEM[local_ind]);
	if (MethodList.OSLMLEM) {
		if (MethodList.MRP) {
			ax[1] += (local_ele * d_MRP_MLEM[local_ind]);
		}
		if (MethodList.Quad) {
			ax[2] += (local_ele * d_Quad_MLEM[local_ind]);
		}
		if (MethodList.L) {
			ax[3] += (local_ele * d_L_MLEM[local_ind]);
		}
		if (MethodList.FMH) {
			ax[4] += (local_ele * d_FMH_MLEM[local_ind]);
		}
		if (MethodList.WeightedMean) {
			ax[5] += (local_ele * d_WeightedMean_MLEM[local_ind]);
		}
		if (MethodList.TV) {
			ax[6] += (local_ele * d_TV_MLEM[local_ind]);
		}
		if (MethodList.AD) {
			ax[7] += (local_ele * d_AD_MLEM[local_ind]);
		}
		if (MethodList.APLS) {
			ax[8] += (local_ele * d_APLS_MLEM[local_ind]);
		}
		if (MethodList.TGV) {
			ax[9] += (local_ele * d_TGV_MLEM[local_ind]);
		}
	}
}

// Nominator (backprojection) in MLEM
inline void nominator_mlem(const RecMethodsOpenCL MethodList, float ax[], const float d_Sino, const float d_epps, const float temp) {
	if (MethodList.MLEM) {
		//nominator_multi(&ax[0], d_Sino, d_epps, temp, randoms_correction, d_sc_ra);
		//if (ax[0] == 0.f)
		//	ax[0] = d_epps;
		//else
		//	ax[0] *= temp;
		//ax[0] = d_Sino / ax[0];
	}
	if (MethodList.OSLMLEM) {
		if (MethodList.MRP) {
			//nominator_multi(&ax[1], d_Sino, d_epps, temp, randoms_correction, d_sc_ra);
			//if (ax[1] == 0.f)
			//	ax[1] = d_epps;
			//else
			//	ax[1] *= temp;
			//ax[1] = d_Sino / ax[1];
		}
		if (MethodList.Quad) {
			//nominator_multi(&ax[2], d_Sino, d_epps, temp, randoms_correction, d_sc_ra);
			//if (ax[2] == 0.f)
			//	ax[2] = d_epps;
			//else
			//	ax[2] *= temp;
			//ax[2] = d_Sino / ax[2];
		}
		if (MethodList.L) {
			//nominator_multi(&ax[3], d_Sino, d_epps, temp, randoms_correction, d_sc_ra);
			//if (ax[3] == 0.f)
			//	ax[3] = d_epps;
			//else
			//	ax[3] *= temp;
			//ax[3] = d_Sino / ax[3];
		}
		if (MethodList.FMH) {
			//nominator_multi(&ax[4], d_Sino, d_epps, temp, randoms_correction, d_sc_ra);
			//if (ax[4] == 0.f)
			//	ax[4] = d_epps;
			//else
			//	ax[4] *= temp;
			//ax[4] = d_Sino / ax[4];
		}
		if (MethodList.WeightedMean) {
			//nominator_multi(&ax[5], d_Sino, d_epps, temp, randoms_correction, d_sc_ra);
			//if (ax[5] == 0.f)
			//	ax[5] = d_epps;
			//else
			//	ax[5] *= temp;
			//ax[5] = d_Sino / ax[5];
		}
		if (MethodList.TV) {
			//nominator_multi(&ax[6], d_Sino, d_epps, temp, randoms_correction, d_sc_ra);
			//if (ax[6] == 0.f)
			//	ax[6] = d_epps;
			//else
			//	ax[6] *= temp;
			//ax[6] = d_Sino / ax[6];
		}
		if (MethodList.AD) {
			//nominator_multi(&ax[7], d_Sino, d_epps, temp, randoms_correction, d_sc_ra);
			//if (ax[7] == 0.f)
			//	ax[7] = d_epps;
			//else
			//	ax[7] *= temp;
			//ax[7] = d_Sino / ax[7];
		}
		if (MethodList.APLS) {
			//nominator_multi(&ax[8], d_Sino, d_epps, temp, randoms_correction, d_sc_ra);
			//if (ax[8] == 0.f)
			//	ax[8] = d_epps;
			//else
			//	ax[8] *= temp;
			//ax[8] = d_Sino / ax[8];
		}
		if (MethodList.TGV) {
			//nominator_multi(&ax[9], d_Sino, d_epps, temp, randoms_correction, d_sc_ra);
			//if (ax[9] == 0.f)
			//	ax[9] = d_epps;
			//else
			//	ax[9] *= temp;
			//ax[9] = d_Sino / ax[9];
		}
	}
}

// Right-hand side of MLEM
inline void rhs_mlem(const RecMethodsOpenCL MethodList, float local_ele, const float ax[], const uint local_ind, __global float* d_rhs_MLEM,
	__global float* d_rhs_MRP_MLEM, __global float* d_rhs_Quad_MLEM, __global float* d_rhs_L_MLEM, __global float* d_rhs_FMH_MLEM, __global float* d_rhs_WeightedMean_MLEM,
	__global float* d_rhs_TV_MLEM, __global float* d_rhs_AD_MLEM, __global float* d_rhs_APLS_MLEM, __global float* d_rhs_TGV_MLEM) {
	if (MethodList.MLEM)
		atomicAdd_g_f(&d_rhs_MLEM[local_ind], (local_ele * ax[0]));
	if (MethodList.OSLMLEM) {
		if (MethodList.MRP)
			atomicAdd_g_f(&d_rhs_MRP_MLEM[local_ind], (local_ele * ax[1]));
		if (MethodList.Quad)
			atomicAdd_g_f(&d_rhs_Quad_MLEM[local_ind], (local_ele * ax[2]));
		if (MethodList.L)
			atomicAdd_g_f(&d_rhs_L_MLEM[local_ind], (local_ele * ax[3]));
		if (MethodList.FMH)
			atomicAdd_g_f(&d_rhs_FMH_MLEM[local_ind], (local_ele * ax[4]));
		if (MethodList.WeightedMean)
			atomicAdd_g_f(&d_rhs_WeightedMean_MLEM[local_ind], (local_ele * ax[5]));
		if (MethodList.TV)
			atomicAdd_g_f(&d_rhs_TV_MLEM[local_ind], (local_ele * ax[6]));
		if (MethodList.AD)
			atomicAdd_g_f(&d_rhs_AD_MLEM[local_ind], (local_ele * ax[7]));
		if (MethodList.APLS)
			atomicAdd_g_f(&d_rhs_APLS_MLEM[local_ind], (local_ele * ax[8]));
		if (MethodList.TGV)
			atomicAdd_g_f(&d_rhs_TGV_MLEM[local_ind], (local_ele * ax[9]));
	}
}

// Denominator (forward projection)
inline void denominator(float local_ele, const RecMethodsOpenCL MethodList, float ax[], const uint local_ind,
	const __global float* d_OSEM, const __global float* d_RAMLA, const __global float* d_MRAMLA, const __global float* d_ROSEM, const __global float* d_RBI, const __global float* d_DRAMA,
	const __global float* d_COSEM, const __global float* d_ACOSEM,
	const __global float* d_MRP_OSEM, const __global float* d_Quad_OSEM, const __global float* d_L_OSEM, const __global float* d_FMH_OSEM, const __global float* d_WeightedMean_OSEM,
	const __global float* d_TV_OSEM, const __global float* d_AD_OSEM, const __global float* d_APLS_OSEM, const __global float* d_TGV_OSEM,
	const __global float* d_MRP_BSREM, const __global float* d_Quad_BSREM, const __global float* d_L_BSREM, const __global float* d_FMH_BSREM, const __global float* d_WeightedMean_BSREM,
	const __global float* d_TV_BSREM, const __global float* d_AD_BSREM, const __global float* d_APLS_BSREM, const __global float* d_TGV_BSREM,
	const __global float* d_MRP_MBSREM, const __global float* d_Quad_MBSREM, const __global float* d_L_MBSREM, const __global float* d_FMH_MBSREM, const __global float* d_WeightedMean_MBSREM,
	const __global float* d_TV_MBSREM, const __global float* d_AD_MBSREM, const __global float* d_APLS_MBSREM, const __global float* d_TGV_MBSREM,
	const __global float* d_MRP_ROSEM, const __global float* d_Quad_ROSEM, const __global float* d_L_ROSEM, const __global float* d_FMH_ROSEM, const __global float* d_WeightedMean_ROSEM,
	const __global float* d_TV_ROSEM, const __global float* d_AD_ROSEM, const __global float* d_APLS_ROSEM, const __global float* d_TGV_ROSEM,
	const __global float* d_MRP_RBI, const __global float* d_Quad_RBI, const __global float* d_L_RBI, const __global float* d_FMH_RBI, const __global float* d_WeightedMean_RBI,
	const __global float* d_TV_RBI, const __global float* d_AD_RBI, const __global float* d_APLS_RBI, const __global float* d_TGV_RBI,
	const __global float* d_MRP_COSEM, const __global float* d_Quad_COSEM, const __global float* d_L_COSEM, const __global float* d_FMH_COSEM, const __global float* d_WeightedMean_COSEM,
	const __global float* d_TV_COSEM, const __global float* d_AD_COSEM, const __global float* d_APLS_COSEM, const __global float* d_TGV_COSEM) {
	//uchar MAP = 0;
	if (MethodList.OSEM == 1 || MethodList.ECOSEM == 1) {
		ax[0] += (local_ele * d_OSEM[local_ind]);
	}
	if (MethodList.MRAMLA == 1)
		ax[1] += (local_ele * d_MRAMLA[local_ind]);
	if (MethodList.RAMLA == 1)
		ax[2] += (local_ele * d_RAMLA[local_ind]);
	if (MethodList.ROSEM == 1)
		ax[3] += (local_ele * d_ROSEM[local_ind]);
	if (MethodList.RBI == 1)
		ax[4] += (local_ele * d_RBI[local_ind]);
	if (MethodList.DRAMA == 1)
		ax[5] += (local_ele * d_DRAMA[local_ind]);
	if (MethodList.COSEM == 1 || MethodList.ECOSEM == 1)
		ax[6] += (local_ele * d_COSEM[local_ind]);
	if (MethodList.ACOSEM == 1)
		ax[7] += (local_ele * d_ACOSEM[local_ind]);
	//if (MAP) {
		if (MethodList.MRP == 1) {
			if (MethodList.OSLOSEM == 1)
				ax[8] += (local_ele * d_MRP_OSEM[local_ind]);
			if (MethodList.BSREM == 1)
				ax[9] += (local_ele * d_MRP_BSREM[local_ind]);
			if (MethodList.MBSREM == 1)
				ax[10] += (local_ele * d_MRP_MBSREM[local_ind]);
			if (MethodList.ROSEMMAP == 1)
				ax[11] += (local_ele * d_MRP_ROSEM[local_ind]);
			if (MethodList.RBIMAP == 1)
				ax[12] += (local_ele * d_MRP_RBI[local_ind]);
			if (MethodList.OSLCOSEM > 0)
				ax[13] += (local_ele * d_MRP_COSEM[local_ind]);
		}
		if (MethodList.Quad == 1) {
			if (MethodList.OSLOSEM == 1)
				ax[14] += (local_ele * d_Quad_OSEM[local_ind]);
			if (MethodList.BSREM == 1)
				ax[15] += (local_ele * d_Quad_BSREM[local_ind]);
			if (MethodList.MBSREM == 1)
				ax[16] += (local_ele * d_Quad_MBSREM[local_ind]);
			if (MethodList.ROSEMMAP == 1)
				ax[17] += (local_ele * d_Quad_ROSEM[local_ind]);
			if (MethodList.RBIMAP == 1)
				ax[18] += (local_ele * d_Quad_RBI[local_ind]);
			if (MethodList.OSLCOSEM > 0)
				ax[19] += (local_ele * d_Quad_COSEM[local_ind]);
		}
		if (MethodList.L == 1) {
			if (MethodList.OSLOSEM == 1)
				ax[20] += (local_ele * d_L_OSEM[local_ind]);
			if (MethodList.BSREM == 1)
				ax[21] += (local_ele * d_L_BSREM[local_ind]);
			if (MethodList.MBSREM == 1)
				ax[22] += (local_ele * d_L_MBSREM[local_ind]);
			if (MethodList.ROSEMMAP == 1)
				ax[23] += (local_ele * d_L_ROSEM[local_ind]);
			if (MethodList.RBIMAP == 1)
				ax[24] += (local_ele * d_L_RBI[local_ind]);
			if (MethodList.OSLCOSEM > 0)
				ax[25] += (local_ele * d_L_COSEM[local_ind]);
		}
		if (MethodList.FMH == 1) {
			if (MethodList.OSLOSEM == 1)
				ax[26] += (local_ele * d_FMH_OSEM[local_ind]);
			if (MethodList.BSREM == 1)
				ax[27] += (local_ele * d_FMH_BSREM[local_ind]);
			if (MethodList.MBSREM == 1)
				ax[28] += (local_ele * d_FMH_MBSREM[local_ind]);
			if (MethodList.ROSEMMAP == 1)
				ax[29] += (local_ele * d_FMH_ROSEM[local_ind]);
			if (MethodList.RBIMAP == 1)
				ax[30] += (local_ele * d_FMH_RBI[local_ind]);
			if (MethodList.OSLCOSEM > 0)
				ax[31] += (local_ele * d_FMH_COSEM[local_ind]);
		}
		if (MethodList.WeightedMean == 1) {
			if (MethodList.OSLOSEM == 1)
				ax[32] += (local_ele * d_WeightedMean_OSEM[local_ind]);
			if (MethodList.BSREM == 1)
				ax[33] += (local_ele * d_WeightedMean_BSREM[local_ind]);
			if (MethodList.MBSREM == 1)
				ax[34] += (local_ele * d_WeightedMean_MBSREM[local_ind]);
			if (MethodList.ROSEMMAP == 1)
				ax[35] += (local_ele * d_WeightedMean_ROSEM[local_ind]);
			if (MethodList.RBIMAP == 1)
				ax[36] += (local_ele * d_WeightedMean_RBI[local_ind]);
			if (MethodList.OSLCOSEM > 0)
				ax[37] += (local_ele * d_WeightedMean_COSEM[local_ind]);
		}
		if (MethodList.TV == 1) {
			if (MethodList.OSLOSEM == 1)
				ax[38] += (local_ele * d_TV_OSEM[local_ind]);
			if (MethodList.BSREM == 1)
				ax[39] += (local_ele * d_TV_BSREM[local_ind]);
			if (MethodList.MBSREM == 1)
				ax[40] += (local_ele * d_TV_MBSREM[local_ind]);
			if (MethodList.ROSEMMAP == 1)
				ax[41] += (local_ele * d_TV_ROSEM[local_ind]);
			if (MethodList.RBIMAP == 1)
				ax[42] += (local_ele * d_TV_RBI[local_ind]);
			if (MethodList.OSLCOSEM > 0)
				ax[43] += (local_ele * d_TV_COSEM[local_ind]);
		}
		if (MethodList.AD == 1) {
			if (MethodList.OSLOSEM == 1)
				ax[44] += (local_ele * d_AD_OSEM[local_ind]);
			if (MethodList.BSREM == 1)
				ax[45] += (local_ele * d_AD_BSREM[local_ind]);
			if (MethodList.MBSREM == 1)
				ax[46] += (local_ele * d_AD_MBSREM[local_ind]);
			if (MethodList.ROSEMMAP == 1)
				ax[47] += (local_ele * d_AD_ROSEM[local_ind]);
			if (MethodList.RBIMAP == 1)
				ax[48] += (local_ele * d_AD_RBI[local_ind]);
			if (MethodList.OSLCOSEM > 0)
				ax[49] += (local_ele * d_AD_COSEM[local_ind]);
		}
		if (MethodList.APLS == 1) {
			if (MethodList.OSLOSEM == 1)
				ax[50] += (local_ele * d_APLS_OSEM[local_ind]);
			if (MethodList.BSREM == 1)
				ax[51] += (local_ele * d_APLS_BSREM[local_ind]);
			if (MethodList.MBSREM == 1)
				ax[52] += (local_ele * d_APLS_MBSREM[local_ind]);
			if (MethodList.ROSEMMAP == 1)
				ax[53] += (local_ele * d_APLS_ROSEM[local_ind]);
			if (MethodList.RBIMAP == 1)
				ax[54] += (local_ele * d_APLS_RBI[local_ind]);
			if (MethodList.OSLCOSEM > 0)
				ax[55] += (local_ele * d_APLS_COSEM[local_ind]);
		}
		if (MethodList.TGV == 1) {
			if (MethodList.OSLOSEM == 1)
				ax[56] += (local_ele * d_TGV_OSEM[local_ind]);
			if (MethodList.BSREM == 1)
				ax[57] += (local_ele * d_TGV_BSREM[local_ind]);
			if (MethodList.MBSREM == 1)
				ax[58] += (local_ele * d_TGV_MBSREM[local_ind]);
			if (MethodList.ROSEMMAP == 1)
				ax[59] += (local_ele * d_TGV_ROSEM[local_ind]);
			if (MethodList.RBIMAP == 1)
				ax[60] += (local_ele * d_TGV_RBI[local_ind]);
			if (MethodList.OSLCOSEM > 0)
				ax[61] += (local_ele * d_TGV_COSEM[local_ind]);
		}
	//}
}

// Nominator (backprojection)
inline void nominator(const RecMethodsOpenCL MethodList, float ax[], const float d_Sino, __constant float* d_epsilon_mramla, const float d_epps, const float temp) {
	if (MethodList.OSEM || MethodList.ECOSEM) {
		if (ax[0] == 0.f)
			ax[0] = d_epps;
		else
			ax[0] *= temp;
		ax[0] = d_Sino / ax[0];
	}
	if (MethodList.MRAMLA) {
		if (ax[1] <= 0.f)
			ax[1] = d_epps;
		else
			ax[1] *= temp;
		if (ax[1] > *d_epsilon_mramla)
			ax[1] = d_Sino / ax[1];
		else
			ax[1] = d_Sino / ax[1] - (d_Sino / native_powr(ax[1], 2)) * (ax[1] - *d_epsilon_mramla);
	}
	if (MethodList.RAMLA) {
		if (ax[2] == 0.f)
			ax[2] = d_epps;
		else
			ax[2] *= temp;
		ax[2] = d_Sino / ax[2];
	}
	if (MethodList.ROSEM) {
		if (ax[3] == 0.f)
			ax[3] = d_epps;
		else
			ax[3] *= temp;
		ax[3] = d_Sino / ax[3];
	}
	if (MethodList.RBI) {
		if (ax[4] == 0.f)
			ax[4] = d_epps;
		else
			ax[4] *= temp;
		ax[4] = d_Sino / ax[4];
	}
	if (MethodList.DRAMA) {
		if (ax[5] == 0.f)
			ax[5] = d_epps;
		else
			ax[5] *= temp;
		ax[5] = d_Sino / ax[5];
	}
	if (MethodList.COSEM || MethodList.ECOSEM) {
		if (ax[6] == 0.f)
			ax[6] = d_epps;
		else
			ax[6] *= temp;
		ax[6] = d_Sino / ax[6];
	}
	if (MethodList.ACOSEM) {
		if (ax[7] == 0.f)
			ax[7] = d_epps;
		else
			ax[7] *= temp;
		ax[7] = d_Sino / ax[7];
	}
	if (MethodList.MRP) {
		if (MethodList.OSLOSEM) {
			if (ax[8] == 0.f)
				ax[8] = d_epps;
			else
				ax[8] *= temp;
			ax[8] = d_Sino / ax[8];
		}
		if (MethodList.BSREM) {
			if (ax[9] == 0.f)
				ax[9] = d_epps;
			else
				ax[9] *= temp;
			ax[9] = d_Sino / ax[9];
		}
		if (MethodList.MBSREM) {
			if (ax[10] == 0.f)
				ax[10] = d_epps;
			else
				ax[10] *= temp;
			ax[10] = d_Sino / ax[10];
		}
		if (MethodList.ROSEMMAP) {
			if (ax[11] == 0.f)
				ax[11] = d_epps;
			else
				ax[11] *= temp;
			ax[11] = d_Sino / ax[11];
		}
		if (MethodList.RBIMAP) {
			if (ax[12] == 0.f)
				ax[12] = d_epps;
			else
				ax[12] *= temp;
			ax[12] = d_Sino / ax[12];
		}
		if (MethodList.OSLCOSEM > 0) {
			if (ax[13] == 0.f)
				ax[13] = d_epps;
			else
				ax[13] *= temp;
			ax[13] = d_Sino / ax[13];
		}
	}
	if (MethodList.Quad) {
		if (MethodList.OSLOSEM) {
			if (ax[14] == 0.f)
				ax[14] = d_epps;
			else
				ax[14] *= temp;
			ax[14] = d_Sino / ax[14];
		}
		if (MethodList.BSREM) {
			if (ax[15] == 0.f)
				ax[15] = d_epps;
			else
				ax[15] *= temp;
			ax[15] = d_Sino / ax[15];
		}
		if (MethodList.MBSREM) {
			if (ax[16] == 0.f)
				ax[16] = d_epps;
			else
				ax[16] *= temp;
			ax[16] = d_Sino / ax[16];
		}
		if (MethodList.ROSEMMAP) {
			if (ax[17] == 0.f)
				ax[17] = d_epps;
			else
				ax[17] *= temp;
			ax[17] = d_Sino / ax[17];
		}
		if (MethodList.RBIMAP) {
			if (ax[18] == 0.f)
				ax[18] = d_epps;
			else
				ax[18] *= temp;
			ax[18] = d_Sino / ax[18];
		}
		if (MethodList.OSLCOSEM > 0) {
			if (ax[19] == 0.f)
				ax[19] = d_epps;
			else
				ax[19] *= temp;
			ax[19] = d_Sino / ax[19];
		}
	}
	if (MethodList.L) {
		if (MethodList.OSLOSEM) {
			if (ax[20] == 0.f)
				ax[20] = d_epps;
			else
				ax[20] *= temp;
			ax[20] = d_Sino / ax[20];
		}
		if (MethodList.BSREM) {
			if (ax[21] == 0.f)
				ax[21] = d_epps;
			else
				ax[21] *= temp;
			ax[21] = d_Sino / ax[21];
		}
		if (MethodList.MBSREM) {
			if (ax[22] == 0.f)
				ax[22] = d_epps;
			else
				ax[22] *= temp;
			ax[22] = d_Sino / ax[22];
		}
		if (MethodList.ROSEMMAP) {
			if (ax[23] == 0.f)
				ax[23] = d_epps;
			else
				ax[23] *= temp;
			ax[23] = d_Sino / ax[23];
		}
		if (MethodList.RBIMAP) {
			if (ax[24] == 0.f)
				ax[24] = d_epps;
			else
				ax[24] *= temp;
			ax[24] = d_Sino / ax[24];
		}
		if (MethodList.OSLCOSEM > 0) {
			if (ax[25] == 0.f)
				ax[25] = d_epps;
			else
				ax[25] *= temp;
			ax[25] = d_Sino / ax[25];
		}
	}
	if (MethodList.FMH) {
		if (MethodList.OSLOSEM) {
			if (ax[26] == 0.f)
				ax[26] = d_epps;
			else
				ax[26] *= temp;
			ax[26] = d_Sino / ax[26];
		}
		if (MethodList.BSREM) {
			if (ax[27] == 0.f)
				ax[27] = d_epps;
			else
				ax[27] *= temp;
			ax[27] = d_Sino / ax[27];
		}
		if (MethodList.MBSREM) {
			if (ax[28] == 0.f)
				ax[28] = d_epps;
			else
				ax[28] *= temp;
			ax[28] = d_Sino / ax[28];
		}
		if (MethodList.ROSEMMAP) {
			if (ax[29] == 0.f)
				ax[29] = d_epps;
			else
				ax[29] *= temp;
			ax[29] = d_Sino / ax[29];
		}
		if (MethodList.RBIMAP) {
			if (ax[30] == 0.f)
				ax[30] = d_epps;
			else
				ax[30] *= temp;
			ax[30] = d_Sino / ax[30];
		}
		if (MethodList.OSLCOSEM > 0) {
			if (ax[31] == 0.f)
				ax[31] = d_epps;
			else
				ax[31] *= temp;
			ax[31] = d_Sino / ax[31];
		}
	}
	if (MethodList.WeightedMean) {
		if (MethodList.OSLOSEM) {
			if (ax[32] == 0.f)
				ax[32] = d_epps;
			else
				ax[32] *= temp;
			ax[32] = d_Sino / ax[32];
		}
		if (MethodList.BSREM) {
			if (ax[33] == 0.f)
				ax[33] = d_epps;
			else
				ax[33] *= temp;
			ax[33] = d_Sino / ax[33];
		}
		if (MethodList.MBSREM) {
			if (ax[34] == 0.f)
				ax[34] = d_epps;
			else
				ax[34] *= temp;
			ax[34] = d_Sino / ax[34];
		}
		if (MethodList.ROSEMMAP) {
			if (ax[35] == 0.f)
				ax[35] = d_epps;
			else
				ax[35] *= temp;
			ax[35] = d_Sino / ax[35];
		}
		if (MethodList.RBIMAP) {
			if (ax[36] == 0.f)
				ax[36] = d_epps;
			else
				ax[36] *= temp;
			ax[36] = d_Sino / ax[36];
		}
		if (MethodList.OSLCOSEM > 0) {
			if (ax[37] == 0.f)
				ax[37] = d_epps;
			else
				ax[37] *= temp;
			ax[37] = d_Sino / ax[37];
		}
	}
	if (MethodList.TV) {
		if (MethodList.OSLOSEM) {
			if (ax[38] == 0.f)
				ax[38] = d_epps;
			else
				ax[38] *= temp;
			ax[38] = d_Sino / ax[38];
		}
		if (MethodList.BSREM) {
			if (ax[39] == 0.f)
				ax[39] = d_epps;
			else
				ax[39] *= temp;
			ax[39] = d_Sino / ax[39];
		}
		if (MethodList.MBSREM) {
			if (ax[40] == 0.f)
				ax[40] = d_epps;
			else
				ax[40] *= temp;
			ax[40] = d_Sino / ax[40];
		}
		if (MethodList.ROSEMMAP) {
			if (ax[41] == 0.f)
				ax[41] = d_epps;
			else
				ax[41] *= temp;
			ax[41] = d_Sino / ax[41];
		}
		if (MethodList.RBIMAP) {
			if (ax[42] == 0.f)
				ax[42] = d_epps;
			else
				ax[42] *= temp;
			ax[42] = d_Sino / ax[42];
		}
		if (MethodList.OSLCOSEM > 0) {
			if (ax[43] == 0.f)
				ax[43] = d_epps;
			else
				ax[43] *= temp;
			ax[43] = d_Sino / ax[43];
		}
	}
	if (MethodList.AD) {
		if (MethodList.OSLOSEM) {
			if (ax[44] == 0.f)
				ax[44] = d_epps;
			else
				ax[44] *= temp;
			ax[44] = d_Sino / ax[44];
		}
		if (MethodList.BSREM) {
			if (ax[45] == 0.f)
				ax[45] = d_epps;
			else
				ax[45] *= temp;
			ax[45] = d_Sino / ax[45];
		}
		if (MethodList.MBSREM) {
			if (ax[46] == 0.f)
				ax[46] = d_epps;
			else
				ax[46] *= temp;
			ax[46] = d_Sino / ax[46];
		}
		if (MethodList.ROSEMMAP) {
			if (ax[47] == 0.f)
				ax[47] = d_epps;
			else
				ax[47] *= temp;
			ax[47] = d_Sino / ax[47];
		}
		if (MethodList.RBIMAP) {
			if (ax[48] == 0.f)
				ax[48] = d_epps;
			else
				ax[48] *= temp;
			ax[48] = d_Sino / ax[48];
		}
		if (MethodList.OSLCOSEM > 0) {
			if (ax[49] == 0.f)
				ax[49] = d_epps;
			else
				ax[49] *= temp;
			ax[49] = d_Sino / ax[49];
		}
	}
	if (MethodList.APLS) {
		if (MethodList.OSLOSEM) {
			if (ax[50] == 0.f)
				ax[50] = d_epps;
			else
				ax[50] *= temp;
			ax[50] = d_Sino / ax[50];
		}
		if (MethodList.BSREM) {
			if (ax[51] == 0.f)
				ax[51] = d_epps;
			else
				ax[51] *= temp;
			ax[51] = d_Sino / ax[51];
		}
		if (MethodList.MBSREM) {
			if (ax[52] == 0.f)
				ax[52] = d_epps;
			else
				ax[52] *= temp;
			ax[52] = d_Sino / ax[52];
		}
		if (MethodList.ROSEMMAP) {
			if (ax[53] == 0.f)
				ax[53] = d_epps;
			else
				ax[53] *= temp;
			ax[53] = d_Sino / ax[53];
		}
		if (MethodList.RBIMAP) {
			if (ax[54] == 0.f)
				ax[54] = d_epps;
			else
				ax[54] *= temp;
			ax[54] = d_Sino / ax[54];
		}
		if (MethodList.OSLCOSEM > 0) {
			if (ax[55] == 0.f)
				ax[55] = d_epps;
			else
				ax[55] *= temp;
			ax[55] = d_Sino / ax[55];
		}
	}
	if (MethodList.TGV) {
		if (MethodList.OSLOSEM) {
			if (ax[56] == 0.f)
				ax[56] = d_epps;
			else
				ax[56] *= temp;
			ax[56] = d_Sino / ax[56];
		}
		if (MethodList.BSREM) {
			if (ax[57] == 0.f)
				ax[57] = d_epps;
			else
				ax[57] *= temp;
			ax[57] = d_Sino / ax[57];
		}
		if (MethodList.MBSREM) {
			if (ax[58] == 0.f)
				ax[58] = d_epps;
			else
				ax[58] *= temp;
			ax[58] = d_Sino / ax[58];
		}
		if (MethodList.ROSEMMAP) {
			if (ax[59] == 0.f)
				ax[59] = d_epps;
			else
				ax[59] *= temp;
			ax[59] = d_Sino / ax[59];
		}
		if (MethodList.RBIMAP) {
			if (ax[60] == 0.f)
				ax[60] = d_epps;
			else
				ax[60] *= temp;
			ax[60] = d_Sino / ax[60];
		}
		if (MethodList.OSLCOSEM > 0) {
			if (ax[61] == 0.f)
				ax[61] = d_epps;
			else
				ax[61] *= temp;
			ax[61] = d_Sino / ax[61];
		}
	}
}

// Right-hand side
inline void rhs(const RecMethodsOpenCL MethodList, float local_ele, const float ax[], const uint local_ind,
	__global float* d_rhs_OSEM, __global float* d_rhs_RAMLA, __global float* d_rhs_MRAMLA, __global float* d_rhs_ROSEM, __global float* d_rhs_RBI, __global float* d_rhs_DRAMA,
	__global float* d_rhs_COSEM, __global float* d_rhs_ACOSEM,
	__global float* d_rhs_MRP_OSEM, __global float* d_rhs_Quad_OSEM, __global float* d_rhs_L_OSEM, __global float* d_rhs_FMH_OSEM, __global float* d_rhs_WeightedMean_OSEM,
	__global float* d_rhs_TV_OSEM, __global float* d_rhs_AD_OSEM, __global float* d_rhs_APLS_OSEM, __global float* d_rhs_TGV_OSEM,
	__global float* d_rhs_MRP_BSREM, __global float* d_rhs_Quad_BSREM, __global float* d_rhs_L_BSREM, __global float* d_rhs_FMH_BSREM, __global float* d_rhs_WeightedMean_BSREM,
	__global float* d_rhs_TV_BSREM, __global float* d_rhs_AD_BSREM, __global float* d_rhs_APLS_BSREM, __global float* d_rhs_TGV_BSREM,
	__global float* d_rhs_MRP_MBSREM, __global float* d_rhs_Quad_MBSREM, __global float* d_rhs_L_MBSREM, __global float* d_rhs_FMH_MBSREM, __global float* d_rhs_WeightedMean_MBSREM,
	__global float* d_rhs_TV_MBSREM, __global float* d_rhs_AD_MBSREM, __global float* d_rhs_APLS_MBSREM, __global float* d_rhs_TGV_MBSREM,
	__global float* d_rhs_MRP_ROSEM, __global float* d_rhs_Quad_ROSEM, __global float* d_rhs_L_ROSEM, __global float* d_rhs_FMH_ROSEM, __global float* d_rhs_WeightedMean_ROSEM,
	__global float* d_rhs_TV_ROSEM, __global float* d_rhs_AD_ROSEM, __global float* d_rhs_APLS_ROSEM, __global float* d_rhs_TGV_ROSEM,
	__global float* d_rhs_MRP_RBI, __global float* d_rhs_Quad_RBI, __global float* d_rhs_L_RBI, __global float* d_rhs_FMH_RBI, __global float* d_rhs_WeightedMean_RBI,
	__global float* d_rhs_TV_RBI, __global float* d_rhs_AD_RBI, __global float* d_rhs_APLS_RBI, __global float* d_rhs_TGV_RBI,
	__global float* d_rhs_MRP_COSEM, __global float* d_rhs_Quad_COSEM, __global float* d_rhs_L_COSEM, __global float* d_rhs_FMH_COSEM, __global float* d_rhs_WeightedMean_COSEM,
	__global float* d_rhs_TV_COSEM, __global float* d_rhs_AD_COSEM, __global float* d_rhs_APLS_COSEM, __global float* d_rhs_TGV_COSEM,
	const float d_h,
	const __global float* d_MRP_COSEM, const __global float* d_Quad_COSEM, const __global float* d_L_COSEM, const __global float* d_FMH_COSEM, const __global float* d_WeightedMean_COSEM,
	const __global float* d_TV_COSEM, const __global float* d_AD_COSEM, const __global float* d_APLS_COSEM, const __global float* d_TGV_COSEM, const __global float* d_COSEM,
	const __global float* d_ACOSEM) {
	if (MethodList.OSEM || MethodList.ECOSEM)
		atomicAdd_g_f(&d_rhs_OSEM[local_ind], (local_ele * ax[0]));
	if (MethodList.MRAMLA)
		atomicAdd_g_f(&d_rhs_MRAMLA[local_ind], (local_ele * ax[1]));
	if (MethodList.RAMLA)
		atomicAdd_g_f(&d_rhs_RAMLA[local_ind], (local_ele * ax[2]));
	if (MethodList.ROSEM)
		atomicAdd_g_f(&d_rhs_ROSEM[local_ind], (local_ele * ax[3]));
	if (MethodList.RBI)
		atomicAdd_g_f(&d_rhs_RBI[local_ind], (local_ele * ax[4]));
	if (MethodList.DRAMA)
		atomicAdd_g_f(&d_rhs_DRAMA[local_ind], (local_ele * ax[5]));
	if (MethodList.COSEM || MethodList.ECOSEM)
		atomicAdd_g_f(&d_rhs_COSEM[local_ind], ax[6] * (d_COSEM[local_ind] * local_ele));
	if (MethodList.ACOSEM)
		atomicAdd_g_f(&d_rhs_ACOSEM[local_ind], ax[7] * (local_ele * native_powr(d_ACOSEM[local_ind], d_h)));
	if (MethodList.MRP) {
		if (MethodList.OSLOSEM)
			atomicAdd_g_f(&d_rhs_MRP_OSEM[local_ind], (local_ele * ax[8]));
		if (MethodList.BSREM)
			atomicAdd_g_f(&d_rhs_MRP_BSREM[local_ind], (local_ele * ax[9]));
		if (MethodList.MBSREM)
			atomicAdd_g_f(&d_rhs_MRP_MBSREM[local_ind], (local_ele * ax[10]));
		if (MethodList.ROSEMMAP)
			atomicAdd_g_f(&d_rhs_MRP_ROSEM[local_ind], (local_ele * ax[11]));
		if (MethodList.RBIMAP)
			atomicAdd_g_f(&d_rhs_MRP_RBI[local_ind], (local_ele * ax[12]));
		if (MethodList.OSLCOSEM == 1)
			atomicAdd_g_f(&d_rhs_MRP_COSEM[local_ind], ax[13] * (local_ele * native_powr(d_MRP_COSEM[local_ind], d_h)));
		else if (MethodList.OSLCOSEM == 2)
			atomicAdd_g_f(&d_rhs_MRP_COSEM[local_ind], ax[13] * (local_ele * d_MRP_COSEM[local_ind]));
	}
	if (MethodList.Quad) {
		if (MethodList.OSLOSEM)
			atomicAdd_g_f(&d_rhs_Quad_OSEM[local_ind], (local_ele * ax[14]));
		if (MethodList.BSREM)
			atomicAdd_g_f(&d_rhs_Quad_BSREM[local_ind], (local_ele * ax[15]));
		if (MethodList.MBSREM)
			atomicAdd_g_f(&d_rhs_Quad_MBSREM[local_ind], (local_ele * ax[16]));
		if (MethodList.ROSEMMAP)
			atomicAdd_g_f(&d_rhs_Quad_ROSEM[local_ind], (local_ele * ax[17]));
		if (MethodList.RBIMAP)
			atomicAdd_g_f(&d_rhs_Quad_RBI[local_ind], (local_ele * ax[18]));
		if (MethodList.OSLCOSEM == 1)
			atomicAdd_g_f(&d_rhs_Quad_COSEM[local_ind], ax[19] * (local_ele * native_powr(d_Quad_COSEM[local_ind], d_h)));
		else if (MethodList.OSLCOSEM == 2)
			atomicAdd_g_f(&d_rhs_Quad_COSEM[local_ind], ax[19] * (local_ele * d_Quad_COSEM[local_ind]));
	}
	if (MethodList.L) {
		if (MethodList.OSLOSEM)
			atomicAdd_g_f(&d_rhs_L_OSEM[local_ind], (local_ele * ax[20]));
		if (MethodList.BSREM)
			atomicAdd_g_f(&d_rhs_L_BSREM[local_ind], (local_ele * ax[21]));
		if (MethodList.MBSREM)
			atomicAdd_g_f(&d_rhs_L_MBSREM[local_ind], (local_ele * ax[22]));
		if (MethodList.ROSEMMAP)
			atomicAdd_g_f(&d_rhs_L_ROSEM[local_ind], (local_ele * ax[23]));
		if (MethodList.RBIMAP)
			atomicAdd_g_f(&d_rhs_L_RBI[local_ind], (local_ele * ax[24]));
		if (MethodList.OSLCOSEM == 1)
			atomicAdd_g_f(&d_rhs_L_COSEM[local_ind], ax[25] * (local_ele * native_powr(d_L_COSEM[local_ind], d_h)));
		else if (MethodList.OSLCOSEM == 2)
			atomicAdd_g_f(&d_rhs_L_COSEM[local_ind], ax[25] * (local_ele * d_L_COSEM[local_ind]));
	}
	if (MethodList.FMH) {
		if (MethodList.OSLOSEM)
			atomicAdd_g_f(&d_rhs_FMH_OSEM[local_ind], (local_ele * ax[26]));
		if (MethodList.BSREM)
			atomicAdd_g_f(&d_rhs_FMH_BSREM[local_ind], (local_ele * ax[27]));
		if (MethodList.MBSREM)
			atomicAdd_g_f(&d_rhs_FMH_MBSREM[local_ind], (local_ele * ax[28]));
		if (MethodList.ROSEMMAP)
			atomicAdd_g_f(&d_rhs_FMH_ROSEM[local_ind], (local_ele * ax[29]));
		if (MethodList.RBIMAP)
			atomicAdd_g_f(&d_rhs_FMH_RBI[local_ind], (local_ele * ax[30]));
		if (MethodList.OSLCOSEM == 1)
			atomicAdd_g_f(&d_rhs_FMH_COSEM[local_ind], ax[31] * (local_ele * native_powr(d_FMH_COSEM[local_ind], d_h)));
		else if (MethodList.OSLCOSEM == 2)
			atomicAdd_g_f(&d_rhs_FMH_COSEM[local_ind], ax[31] * (local_ele * d_FMH_COSEM[local_ind]));
	}
	if (MethodList.WeightedMean) {
		if (MethodList.OSLOSEM)
			atomicAdd_g_f(&d_rhs_WeightedMean_OSEM[local_ind], (local_ele * ax[32]));
		if (MethodList.BSREM)
			atomicAdd_g_f(&d_rhs_WeightedMean_BSREM[local_ind], (local_ele * ax[33]));
		if (MethodList.MBSREM)
			atomicAdd_g_f(&d_rhs_WeightedMean_MBSREM[local_ind], (local_ele * ax[34]));
		if (MethodList.ROSEMMAP)
			atomicAdd_g_f(&d_rhs_WeightedMean_ROSEM[local_ind], (local_ele * ax[35]));
		if (MethodList.RBIMAP)
			atomicAdd_g_f(&d_rhs_WeightedMean_RBI[local_ind], (local_ele * ax[36]));
		if (MethodList.OSLCOSEM == 1)
			atomicAdd_g_f(&d_rhs_WeightedMean_COSEM[local_ind], ax[37] * (local_ele * native_powr(d_WeightedMean_COSEM[local_ind], d_h)));
		else if (MethodList.OSLCOSEM == 2)
			atomicAdd_g_f(&d_rhs_WeightedMean_COSEM[local_ind], ax[37] * (local_ele * d_WeightedMean_COSEM[local_ind]));
	}
	if (MethodList.TV) {
		if (MethodList.OSLOSEM)
			atomicAdd_g_f(&d_rhs_TV_OSEM[local_ind], (local_ele * ax[38]));
		if (MethodList.BSREM)
			atomicAdd_g_f(&d_rhs_TV_BSREM[local_ind], (local_ele * ax[39]));
		if (MethodList.MBSREM)
			atomicAdd_g_f(&d_rhs_TV_MBSREM[local_ind], (local_ele * ax[40]));
		if (MethodList.ROSEMMAP)
			atomicAdd_g_f(&d_rhs_TV_ROSEM[local_ind], (local_ele * ax[41]));
		if (MethodList.RBIMAP)
			atomicAdd_g_f(&d_rhs_TV_RBI[local_ind], (local_ele * ax[42]));
		if (MethodList.OSLCOSEM == 1)
			atomicAdd_g_f(&d_rhs_TV_COSEM[local_ind], ax[43] * (local_ele * native_powr(d_TV_COSEM[local_ind], d_h)));
		else if (MethodList.OSLCOSEM == 2)
			atomicAdd_g_f(&d_rhs_TV_COSEM[local_ind], ax[43] * (local_ele * d_TV_COSEM[local_ind]));
	}
	if (MethodList.AD) {
		if (MethodList.OSLOSEM)
			atomicAdd_g_f(&d_rhs_AD_OSEM[local_ind], (local_ele * ax[44]));
		if (MethodList.BSREM)
			atomicAdd_g_f(&d_rhs_AD_BSREM[local_ind], (local_ele * ax[45]));
		if (MethodList.MBSREM)
			atomicAdd_g_f(&d_rhs_AD_MBSREM[local_ind], (local_ele * ax[46]));
		if (MethodList.ROSEMMAP)
			atomicAdd_g_f(&d_rhs_AD_ROSEM[local_ind], (local_ele * ax[47]));
		if (MethodList.RBIMAP)
			atomicAdd_g_f(&d_rhs_AD_RBI[local_ind], (local_ele * ax[48]));
		if (MethodList.OSLCOSEM == 1)
			atomicAdd_g_f(&d_rhs_AD_COSEM[local_ind], ax[49] * (local_ele * native_powr(d_AD_COSEM[local_ind], d_h)));
		else if (MethodList.OSLCOSEM == 2)
			atomicAdd_g_f(&d_rhs_AD_COSEM[local_ind], ax[49] * (local_ele * d_AD_COSEM[local_ind]));
	}
	if (MethodList.APLS) {
		if (MethodList.OSLOSEM)
			atomicAdd_g_f(&d_rhs_APLS_OSEM[local_ind], (local_ele * ax[50]));
		if (MethodList.BSREM)
			atomicAdd_g_f(&d_rhs_APLS_BSREM[local_ind], (local_ele * ax[51]));
		if (MethodList.MBSREM)
			atomicAdd_g_f(&d_rhs_APLS_MBSREM[local_ind], (local_ele * ax[52]));
		if (MethodList.ROSEMMAP)
			atomicAdd_g_f(&d_rhs_APLS_ROSEM[local_ind], (local_ele * ax[53]));
		if (MethodList.RBIMAP)
			atomicAdd_g_f(&d_rhs_APLS_RBI[local_ind], (local_ele * ax[54]));
		if (MethodList.OSLCOSEM == 1)
			atomicAdd_g_f(&d_rhs_APLS_COSEM[local_ind], ax[55] * (local_ele * native_powr(d_APLS_COSEM[local_ind], d_h)));
		else if (MethodList.OSLCOSEM == 2)
			atomicAdd_g_f(&d_rhs_APLS_COSEM[local_ind], ax[55] * (local_ele * d_APLS_COSEM[local_ind]));
	}
	if (MethodList.TGV) {
		if (MethodList.OSLOSEM)
			atomicAdd_g_f(&d_rhs_TGV_OSEM[local_ind], (local_ele * ax[56]));
		if (MethodList.BSREM)
			atomicAdd_g_f(&d_rhs_TGV_BSREM[local_ind], (local_ele * ax[57]));
		if (MethodList.MBSREM)
			atomicAdd_g_f(&d_rhs_TGV_MBSREM[local_ind], (local_ele * ax[58]));
		if (MethodList.ROSEMMAP)
			atomicAdd_g_f(&d_rhs_TGV_ROSEM[local_ind], (local_ele * ax[59]));
		if (MethodList.RBIMAP)
			atomicAdd_g_f(&d_rhs_TGV_RBI[local_ind], (local_ele * ax[60]));
		if (MethodList.OSLCOSEM == 1)
			atomicAdd_g_f(&d_rhs_TGV_COSEM[local_ind], ax[61] * (local_ele * native_powr(d_TGV_COSEM[local_ind], d_h)));
		else if (MethodList.OSLCOSEM == 2)
			atomicAdd_g_f(&d_rhs_TGV_COSEM[local_ind], ax[61] * (local_ele * d_TGV_COSEM[local_ind]));
	}
}

// Get the detector coordinates for the current (raw) measurement
inline void get_detector_coordinates_raw(const __global float *d_x, const __global float *d_y, const __global float *d_zdet, const __global ushort* d_L, const uint d_det_per_ring,
	const uint idx, __constant uint *d_pseudos, const uint d_pRows, float *xs, float* xd, float* ys, float* yd, float* zs, float* zd) {
	uint ps = 0;
	// Get the current detector numbers
	const uint detektorit1 = convert_uint(d_L[idx * 2u]) - 1u;
	const uint detektorit2 = convert_uint(d_L[idx * 2u + 1u]) - 1u;
	// Which ring
	const uint loop1 = ((detektorit1) / d_det_per_ring);
	const uint loop2 = ((detektorit2) / d_det_per_ring);
	// same ring
	if (loop1 == loop2) {
		*zs = d_zdet[loop1];
		*zd = *zs;
	}
	else {
		*zs = d_zdet[loop1];
		*zd = d_zdet[loop2];
	}
	// if a pseudo ring is present, move the z-detector coordinate beyond the pseudo ring (z_det includes also pseudo coordinates)
	if (loop1 >= d_pseudos[0]) {
		// loop through all the pseudo rings
		for (uint kk = 0u; kk < d_pRows; kk++) {
			ps = 1u;
			// First pseudo rings
			if (kk + 1u < d_pRows) {
				if (loop1 >= d_pseudos[kk] && loop1 < d_pseudos[kk + 1u]) {
					*zs = d_zdet[loop1 + ps];
					break;
				}
				// move to next
				else
					ps++;
			}
			else {
				// Last pseudo ring passed
				if (loop1 >= d_pseudos[kk])
					*zs = d_zdet[loop1 + ps];
			}
		}
	}
	if (loop2 >= d_pseudos[0]) {
		for (uint kk = 0u; kk < d_pRows; kk++) {
			ps = 1u;
			if (kk + 1u < d_pRows) {
				if (loop2 >= d_pseudos[kk] && loop2 < d_pseudos[kk + 1u]) {
					*zd = d_zdet[loop2 + ps];
					break;
				}
				else
					ps++;
			}
			else {
				if (loop2 >= d_pseudos[kk])
					*zd = d_zdet[loop2 + ps];
			}
		}
	}
	// Get the current x- and y-detector coordinates
	*xs = d_x[detektorit1 - d_det_per_ring * (loop1)];
	*xd = d_x[detektorit2 - d_det_per_ring * (loop2)];
	*ys = d_y[detektorit1 - d_det_per_ring * (loop1)];
	*yd = d_y[detektorit2 - d_det_per_ring * (loop2)];
}

// Get the detector coordinates for the current sinogram bin
inline void get_detector_coordinates(const __global uint *d_xyindex, const __global ushort *d_zindex, const uint d_size_x,
	const uint idx, const ushort d_TotSinos, float *xs, float* xd, float* ys, float* yd, float* zs, float* zd, const __global float *d_x, const __global float *d_y,
	const __global float *d_zdet) {

	if (d_xyindex[idx] >= d_size_x) {
		*xs = d_x[d_xyindex[idx]];
		*xd = d_x[d_xyindex[idx] - d_size_x];
		*ys = d_y[d_xyindex[idx]];
		*yd = d_y[d_xyindex[idx] - d_size_x];
	}
	else {
		*xs = d_x[d_xyindex[idx]];
		*xd = d_x[d_xyindex[idx] + d_size_x];
		*ys = d_y[d_xyindex[idx]];
		*yd = d_y[d_xyindex[idx] + d_size_x];
	}
	if (d_zindex[idx] >= d_TotSinos) {
		*zs = d_zdet[d_zindex[idx]];
		*zd = d_zdet[d_zindex[idx] - d_TotSinos];
	}
	else {
		*zs = d_zdet[d_zindex[idx]];
		*zd = d_zdet[d_zindex[idx] + d_TotSinos];
	}
}

// Get the detector coordinates for the current measurement (precomputation phase)
inline void get_detector_coordinates_precomp(const uint d_size_x, const uint idx, const ushort d_TotSinos, float *xs, float* xd, float* ys, float* yd, float* zs,
	float* zd, const __global float *d_x, const __global float *d_y, const __global float *d_zdet) {

	const uint id = idx % d_size_x;
	const uint idz = idx / d_size_x;
	*xs = d_x[id];
	*xd = d_x[id + d_size_x];
	*ys = d_y[id];
	*yd = d_y[id + d_size_x];
	*zs = d_zdet[idz];
	*zd = d_zdet[idz + d_TotSinos];
}

// Compute the voxel index where the current perpendicular measurement starts
inline int perpendicular_start(const float d_b, const float d, const float d_d, const uint d_N) {
	int tempi;
	float start = d_b - d + d_d;
	for (uint ii = 0u; ii < d_N; ii++) {
		if (start > 0.f) {
			tempi = convert_int(ii);
			break;
		}
		start += d_d;
	}
	return tempi;
}

// Compute the probability for the perpendicular elements
inline void perpendicular_elements(const float d_b, const float d_d1, const float d_N1, const float d, const float d_d2, const float d_N2, const __global float* d_atten,
	float* templ_ijk, uint* tempk, const uint d_attenuation_correction, const uint z_loop, const uint d_N, const uint d_NN, const uint d_normalization, const float d_norm) {
	int apu = perpendicular_start(d_b, d, d_d1, d_N1);
	//float start = d_b - d + d_d1;
	//// Find the closest y-index value by finding the smallest y-distance between detector 2 and all the y-pixel coordinates
	//for (int ii = 0; ii < d_N1; ii++) {
	//	if (start > 0.f) {
	//		apu = ii;
	//		break;
	//	}
	//	start += d_d1;
	//}
	//*templ_ijk = d_d2;
	*tempk = convert_uint_sat(apu) * d_N + z_loop * d_N1 * d_N2;
	float temp = d_d2 * convert_float(d_N2);
	// Probability
	temp = 1.f / temp;
	if (d_attenuation_correction == 1u) {
		float jelppi = 0.f;
		for (uint iii = 0u; iii < d_N2; iii++) {
			jelppi += (*templ_ijk * (-d_atten[*tempk + iii * d_NN]));
		}
		temp *= native_exp(jelppi);
	}
	if (d_normalization == 1u)
		temp *= d_norm;
	*templ_ijk = temp * d_d2;
}

// Compute functions (9) and (29) (detector larger than source)
inline void d_g_s(const float tmin, const float t_min, uint* v_min, float* t_0, int* v_u, const float diff, const float b, const float d, const float s) {

	if (tmin == t_min)
		// (11)
		*v_min = 1u;
	else {
		// (2) and (19)
		const float p_t = s + tmin * (diff);
		// (12)
		*v_min = convert_uint_sat_rtp((p_t - b) / d);
	}
	//if (tmax == t_max)
	//	// (13)
	//	*v_max = N;
	//else {
	//	// (2) and (19)
	//	p_t = s + tmax * (diff);
	//	// (14)
	//	*v_max = convert_int_rtz((p_t - b) / d);
	//}
	// (9)
	*t_0 += ((convert_float(*v_min) * d) / (diff));
	//  (29)
	*v_u = 1;
}

// Compute functions (9) and (29) (source larger than detector)
inline void s_g_d(const float tmin, const float t_min, uint* v_max, float* t_0, int* v_u, const float diff, const float b, const float d, const float s, const uint N) {

	if (tmin == t_min)
		// (15)
		*v_max = N - 1u;
	else {
		// (2) and (19)
		const float p_t = s + tmin * (diff);
		// (16)
		*v_max = convert_uint_sat_rtz((p_t - b) / d);
	}
	//if (tmax == t_max)
	//	// (17)
	//	*v_min = 0;
	//else {
	//	// (2) and (19)
	//	p_t = s + tmax * (diff);
	//	// (18)
	//	*v_min = convert_int_rtp((p_t - b) / d);
	//}
	// (9)
	*t_0 += ((convert_float(*v_max) * d) / (diff));
	// (29)
	*v_u = -1;
}

// same as above, but for precomputation phase
inline void d_g_s_precomp(const float tmin, const float t_min, const float tmax, const float t_max, uint* v_min, uint* v_max, float* t_0, int* v_u, const float diff,
	const float b, const float d, const float s, const uint N) {

	if (tmin == t_min)
		// (11)
		*v_min = 1u;
	else {
		// (2) and (19)
		const float p_t = s + tmin * (diff);
		// (12)
		*v_min = convert_uint_sat_rtp((p_t - b) / d);
	}
	if (tmax == t_max)
		// (13)
		*v_max = N;
	else {
		// (2) and (19)
		const float p_t = s + tmax * (diff);
		// (14)
		*v_max = convert_uint_sat_rtz((p_t - b) / d);
	}
	// (9)
	*t_0 += ((convert_float(*v_min) * d) / (diff));
	//  (29)
	*v_u = 1;
}

// same as above, but for precomputation phase
inline void s_g_d_precomp(const float tmin, const float t_min, const float tmax, const float t_max, uint* v_min, uint* v_max, float* t_0, int* v_u, const float diff,
	const float b, const float d, const float s, const uint N) {

	if (tmin == t_min)
		// (15)
		*v_max = N - 1u;
	else {
		// (2) and (19)
		const float p_t = s + tmin * (diff);
		// (16)
		*v_max = convert_uint_sat_rtz((p_t - b) / d);
	}
	if (tmax == t_max)
		// (17)
		*v_min = 0u;
	else {
		// (2) and (19)
		const float p_t = s + tmax * (diff);
		// (18)
		*v_min = convert_uint_rtp((p_t - b) / d);
	}
	// (9)
	*t_0 += ((convert_float(*v_max) * d) / (diff));
	// (29)
	*v_u = -1;
}

// Compute the index of the current voxel
inline uint compute_ind(const int tempj, const int tempi, const int tempk, const uint d_N1, const uint d_N2, const uint d_N, const uint d_Nx, 
	const uint d_Nyx) {
	uint local_ind = convert_uint_sat(tempj) * d_Nx + convert_uint_sat(tempi) + convert_uint_sat(tempk) * d_Nyx;
	if (local_ind >= d_N) {
		if (local_ind - d_N1 >= d_N)
			local_ind -= (d_N1 * d_N2);
		else if (local_ind - 1u >= d_N)
			local_ind -= d_N1;
		else
			local_ind--;
	}
	//else if (local_ind < 0) {
	//	if (local_ind + d_N1 < 0)
	//		local_ind += (d_N1 * d_N2);
	//	else if (local_ind + 1 < 0)
	//		local_ind += d_N1;
	//	else
	//		local_ind++;
	//}
	return local_ind;
}

// compute the distance that the ray traverses in the current voxel
inline float compute_element(float* t0, float* tc, const float L, const float tu, const int u, int* temp_ijk, float* temp) {
	float local_ele = (*t0 - *tc) * L;
	*temp_ijk += u;
	tc = t0;
	*t0 += tu;
	*temp += local_ele;
	return local_ele;
}

inline float compute_matrix_element(const float t0, const float tc, const float L) {
	return (t0 - tc) * L;
}

inline void compute_attenuation(float* tc, float* jelppi, const float LL, const float t0, const int tempi, const int tempj, const int tempk, const uint Nx, const uint Nyx, 
	const __global float* d_atten) {
	*jelppi += (compute_matrix_element(t0, *tc, LL) * -d_atten[tempi + tempj * Nx + Nyx * tempk]);
	*tc = t0;
}

// compute the probability of emission in the current voxel
inline float compute_element_2nd(float* t0, float* tc, const float L, const float tu, const int u, int* temp_ijk, const float temp) {
	float local_ele = (*t0 - *tc) * L * temp;
	*temp_ijk += u;
	tc = t0;
	*t0 += tu;
	return local_ele;
}

// compute the initial voxel index (beginning of the ray)
inline int voxel_index(const float pt, const float diff, const float d, const float apu) {
	return convert_int_rtz((pt * diff - apu) / d);
}

// Denominator (forward projection), multi-GPU version
inline void denominator_multi(float local_ele, float *axOSEM, const __global float* d_OSEM) {
	*axOSEM += (local_ele * *d_OSEM);
}

// Nominator (backprojection), multi-GPU version
inline void nominator_multi(float *axOSEM, const float d_Sino, const float d_epps, const float temp, const uint randoms_correction, const float d_sc_ra) {
	if (*axOSEM == 0.f)
		*axOSEM = d_epps;
	else
		*axOSEM *= temp;
	if (randoms_correction == 1u)
		*axOSEM += d_sc_ra;
	*axOSEM = d_Sino / *axOSEM;
}

// RHS, multi-GPU version
inline void rhs_multi(float local_ele, const float *yaxOSEM, __global float* d_rhs_OSEM) {
	atomicAdd_g_f(d_rhs_OSEM, (local_ele * *yaxOSEM));
}

// compute voxel index, orthogonal distance based ray tracer
inline uint compute_ind_orth(const int tempi, const uint temp_ijk, const uint d_N1, const uint d_N2, const uint d_N) {
	uint local_ind = convert_uint_sat(tempi) * d_N + temp_ijk;
	//if (local_ind >= d_N) {
	//	if (local_ind - d_N1 >= d_N)
	//		local_ind -= (d_N1 * d_N2);
	//	else
	//		local_ind--;
	//}
	//else if (local_ind < 0) {
	//	if (local_ind + d_N1 < 0)
	//		local_ind += (d_N1 * d_N2);
	//	else
	//		local_ind++;
	//}
	return local_ind;
}

// compute orthogonal distance
inline float compute_element_orth(const float x_diff, const float y_diff, const float x_center, const float length_) {
	float local_ele = 1.f - fabs(x_diff + y_diff * x_center) / length_;
	return local_ele;
}

// Denominator (forward projection), orthogonal distance based ray tracer
inline void orth_distance_denominator(const int tempi, const uint Nx, const float y_diff, const float x_diff, __constant float* y_center, __constant float* x_center, const float kerroin,
	const float length_, float* temp, const uint temp_ijk, const RecMethodsOpenCL MethodList, 
	const float local_sino, float ax[], const uint d_Ny, const uint d_N, const int tempj,
	const __global float* d_OSEM, const __global float* d_RAMLA, const __global float* d_MRAMLA, const __global float* d_ROSEM, const __global float* d_RBI, const __global float* d_DRAMA,
	const __global float* d_COSEM, const __global float* d_ACOSEM,
	const __global float* d_MRP_OSEM, const __global float* d_Quad_OSEM, const __global float* d_L_OSEM, const __global float* d_FMH_OSEM, const __global float* d_WeightedMean_OSEM,
	const __global float* d_TV_OSEM, const __global float* d_AD_OSEM, const __global float* d_APLS_OSEM, const __global float* d_TGV_OSEM,
	const __global float* d_MRP_BSREM, const __global float* d_Quad_BSREM, const __global float* d_L_BSREM, const __global float* d_FMH_BSREM, const __global float* d_WeightedMean_BSREM,
	const __global float* d_TV_BSREM, const __global float* d_AD_BSREM, const __global float* d_APLS_BSREM, const __global float* d_TGV_BSREM,
	const __global float* d_MRP_MBSREM, const __global float* d_Quad_MBSREM, const __global float* d_L_MBSREM, const __global float* d_FMH_MBSREM, const __global float* d_WeightedMean_MBSREM,
	const __global float* d_TV_MBSREM, const __global float* d_AD_MBSREM, const __global float* d_APLS_MBSREM, const __global float* d_TGV_MBSREM,
	const __global float* d_MRP_ROSEM, const __global float* d_Quad_ROSEM, const __global float* d_L_ROSEM, const __global float* d_FMH_ROSEM, const __global float* d_WeightedMean_ROSEM,
	const __global float* d_TV_ROSEM, const __global float* d_AD_ROSEM, const __global float* d_APLS_ROSEM, const __global float* d_TGV_ROSEM,
	const __global float* d_MRP_RBI, const __global float* d_Quad_RBI, const __global float* d_L_RBI, const __global float* d_FMH_RBI, const __global float* d_WeightedMean_RBI,
	const __global float* d_TV_RBI, const __global float* d_AD_RBI, const __global float* d_APLS_RBI, const __global float* d_TGV_RBI,
	const __global float* d_MRP_COSEM, const __global float* d_Quad_COSEM, const __global float* d_L_COSEM, const __global float* d_FMH_COSEM, const __global float* d_WeightedMean_COSEM,
	const __global float* d_TV_COSEM, const __global float* d_AD_COSEM, const __global float* d_APLS_COSEM, const __global float* d_TGV_COSEM) {

	const float diff = kerroin - x_diff * y_center[tempj];
	for (int uu = tempi; uu >= 0; uu--) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		*temp += local_ele;
		if (local_sino > 0.f) {
			denominator(local_ele, MethodList, ax, local_ind,
				d_OSEM, d_RAMLA, d_MRAMLA, d_ROSEM, d_RBI, d_DRAMA, d_COSEM, d_ACOSEM,
				d_MRP_OSEM, d_Quad_OSEM, d_L_OSEM, d_FMH_OSEM, d_WeightedMean_OSEM, d_TV_OSEM, d_AD_OSEM, d_APLS_OSEM, d_TGV_OSEM,
				d_MRP_BSREM, d_Quad_BSREM, d_L_BSREM, d_FMH_BSREM, d_WeightedMean_BSREM, d_TV_BSREM, d_AD_BSREM, d_APLS_BSREM, d_TGV_BSREM,
				d_MRP_MBSREM, d_Quad_MBSREM, d_L_MBSREM, d_FMH_MBSREM, d_WeightedMean_MBSREM, d_TV_MBSREM, d_AD_MBSREM, d_APLS_MBSREM, d_TGV_MBSREM,
				d_MRP_ROSEM, d_Quad_ROSEM, d_L_ROSEM, d_FMH_ROSEM, d_WeightedMean_ROSEM, d_TV_ROSEM, d_AD_ROSEM, d_APLS_ROSEM, d_TGV_ROSEM,
				d_MRP_RBI, d_Quad_RBI, d_L_RBI, d_FMH_RBI, d_WeightedMean_RBI, d_TV_RBI, d_AD_RBI, d_APLS_RBI, d_TGV_RBI,
				d_MRP_COSEM, d_Quad_COSEM, d_L_COSEM, d_FMH_COSEM, d_WeightedMean_COSEM, d_TV_COSEM, d_AD_COSEM, d_APLS_COSEM, d_TGV_COSEM);
		}
	}
	for (int uu = tempi + 1; uu < Nx; uu++) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		*temp += local_ele;
		if (local_sino > 0.f) {
			denominator(local_ele, MethodList, ax, local_ind,
				d_OSEM, d_RAMLA, d_MRAMLA, d_ROSEM, d_RBI, d_DRAMA, d_COSEM, d_ACOSEM,
				d_MRP_OSEM, d_Quad_OSEM, d_L_OSEM, d_FMH_OSEM, d_WeightedMean_OSEM, d_TV_OSEM, d_AD_OSEM, d_APLS_OSEM, d_TGV_OSEM,
				d_MRP_BSREM, d_Quad_BSREM, d_L_BSREM, d_FMH_BSREM, d_WeightedMean_BSREM, d_TV_BSREM, d_AD_BSREM, d_APLS_BSREM, d_TGV_BSREM,
				d_MRP_MBSREM, d_Quad_MBSREM, d_L_MBSREM, d_FMH_MBSREM, d_WeightedMean_MBSREM, d_TV_MBSREM, d_AD_MBSREM, d_APLS_MBSREM, d_TGV_MBSREM,
				d_MRP_ROSEM, d_Quad_ROSEM, d_L_ROSEM, d_FMH_ROSEM, d_WeightedMean_ROSEM, d_TV_ROSEM, d_AD_ROSEM, d_APLS_ROSEM, d_TGV_ROSEM,
				d_MRP_RBI, d_Quad_RBI, d_L_RBI, d_FMH_RBI, d_WeightedMean_RBI, d_TV_RBI, d_AD_RBI, d_APLS_RBI, d_TGV_RBI,
				d_MRP_COSEM, d_Quad_COSEM, d_L_COSEM, d_FMH_COSEM, d_WeightedMean_COSEM, d_TV_COSEM, d_AD_COSEM, d_APLS_COSEM, d_TGV_COSEM);
		}
	}
}

// RHS, orthogonal distance based ray tracer
inline void orth_distance_rhs(const int tempi, const uint Nx, const float y_diff, const float x_diff, const float y_center, __constant float* x_center, const float kerroin,
	const float length_, const float temp, const uint temp_ijk, const RecMethodsOpenCL MethodList,
	const float ax[], const uint d_Ny, const uint d_N,
	__global float* d_rhs_OSEM, __global float* d_rhs_RAMLA, __global float* d_rhs_MRAMLA, __global float* d_rhs_ROSEM, __global float* d_rhs_RBI, __global float* d_rhs_DRAMA,
	__global float* d_rhs_COSEM, __global float* d_rhs_ACOSEM,
	__global float* d_rhs_MRP_OSEM, __global float* d_rhs_Quad_OSEM, __global float* d_rhs_L_OSEM, __global float* d_rhs_FMH_OSEM, __global float* d_rhs_WeightedMean_OSEM,
	__global float* d_rhs_TV_OSEM, __global float* d_rhs_AD_OSEM, __global float* d_rhs_APLS_OSEM, __global float* d_rhs_TGV_OSEM,
	__global float* d_rhs_MRP_BSREM, __global float* d_rhs_Quad_BSREM, __global float* d_rhs_L_BSREM, __global float* d_rhs_FMH_BSREM, __global float* d_rhs_WeightedMean_BSREM,
	__global float* d_rhs_TV_BSREM, __global float* d_rhs_AD_BSREM, __global float* d_rhs_APLS_BSREM, __global float* d_rhs_TGV_BSREM,
	__global float* d_rhs_MRP_MBSREM, __global float* d_rhs_Quad_MBSREM, __global float* d_rhs_L_MBSREM, __global float* d_rhs_FMH_MBSREM, __global float* d_rhs_WeightedMean_MBSREM,
	__global float* d_rhs_TV_MBSREM, __global float* d_rhs_AD_MBSREM, __global float* d_rhs_APLS_MBSREM, __global float* d_rhs_TGV_MBSREM,
	__global float* d_rhs_MRP_ROSEM, __global float* d_rhs_Quad_ROSEM, __global float* d_rhs_L_ROSEM, __global float* d_rhs_FMH_ROSEM, __global float* d_rhs_WeightedMean_ROSEM,
	__global float* d_rhs_TV_ROSEM, __global float* d_rhs_AD_ROSEM, __global float* d_rhs_APLS_ROSEM, __global float* d_rhs_TGV_ROSEM,
	__global float* d_rhs_MRP_RBI, __global float* d_rhs_Quad_RBI, __global float* d_rhs_L_RBI, __global float* d_rhs_FMH_RBI, __global float* d_rhs_WeightedMean_RBI,
	__global float* d_rhs_TV_RBI, __global float* d_rhs_AD_RBI, __global float* d_rhs_APLS_RBI, __global float* d_rhs_TGV_RBI,
	__global float* d_rhs_MRP_COSEM, __global float* d_rhs_Quad_COSEM, __global float* d_rhs_L_COSEM, __global float* d_rhs_FMH_COSEM, __global float* d_rhs_WeightedMean_COSEM,
	__global float* d_rhs_TV_COSEM, __global float* d_rhs_AD_COSEM, __global float* d_rhs_APLS_COSEM, __global float* d_rhs_TGV_COSEM,
	const float d_h,
	const __global float* d_MRP_COSEM, const __global float* d_Quad_COSEM, const __global float* d_L_COSEM, const __global float* d_FMH_COSEM, const __global float* d_WeightedMean_COSEM,
	const __global float* d_TV_COSEM, const __global float* d_AD_COSEM, const __global float* d_APLS_COSEM, const __global float* d_TGV_COSEM, const __global float* d_COSEM,
	const __global float* d_ACOSEM, 
	__global float* Summ) {

	const float diff = kerroin - x_diff * y_center;
	for (int uu = tempi; uu >= 0; uu--) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		local_ele *= temp;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		atomicAdd_g_f(&Summ[local_ind], local_ele);
		rhs(MethodList, local_ele, ax, local_ind,
			d_rhs_OSEM, d_rhs_RAMLA, d_rhs_MRAMLA, d_rhs_ROSEM, d_rhs_RBI, d_rhs_DRAMA, d_rhs_COSEM, d_rhs_ACOSEM,
			d_rhs_MRP_OSEM, d_rhs_Quad_OSEM, d_rhs_L_OSEM, d_rhs_FMH_OSEM, d_rhs_WeightedMean_OSEM, d_rhs_TV_OSEM, d_rhs_AD_OSEM, d_rhs_APLS_OSEM, d_rhs_TGV_OSEM,
			d_rhs_MRP_BSREM, d_rhs_Quad_BSREM, d_rhs_L_BSREM, d_rhs_FMH_BSREM, d_rhs_WeightedMean_BSREM, d_rhs_TV_BSREM, d_rhs_AD_BSREM, d_rhs_APLS_BSREM, 
			d_rhs_TGV_BSREM,
			d_rhs_MRP_MBSREM, d_rhs_Quad_MBSREM, d_rhs_L_MBSREM, d_rhs_FMH_MBSREM, d_rhs_WeightedMean_MBSREM, d_rhs_TV_MBSREM, d_rhs_AD_MBSREM, d_rhs_APLS_MBSREM, 
			d_rhs_TGV_MBSREM,
			d_rhs_MRP_ROSEM, d_rhs_Quad_ROSEM, d_rhs_L_ROSEM, d_rhs_FMH_ROSEM, d_rhs_WeightedMean_ROSEM, d_rhs_TV_ROSEM, d_rhs_AD_ROSEM, d_rhs_APLS_ROSEM, 
			d_rhs_TGV_ROSEM,
			d_rhs_MRP_RBI, d_rhs_Quad_RBI, d_rhs_L_RBI, d_rhs_FMH_RBI, d_rhs_WeightedMean_RBI, d_rhs_TV_RBI, d_rhs_AD_RBI, d_rhs_APLS_RBI, d_rhs_TGV_RBI,
			d_rhs_MRP_COSEM, d_rhs_Quad_COSEM, d_rhs_L_COSEM, d_rhs_FMH_COSEM, d_rhs_WeightedMean_COSEM, d_rhs_TV_COSEM, d_rhs_AD_COSEM, d_rhs_APLS_COSEM, 
			d_rhs_TGV_COSEM, d_h,
			d_MRP_COSEM, d_Quad_COSEM, d_L_COSEM, d_FMH_COSEM, d_WeightedMean_COSEM, d_TV_COSEM, d_AD_COSEM, d_APLS_COSEM, d_TGV_COSEM, d_COSEM, d_ACOSEM);
	}
	for (int uu = tempi + 1; uu < Nx; uu++) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		local_ele *= temp;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		rhs(MethodList, local_ele, ax, local_ind,
			d_rhs_OSEM, d_rhs_RAMLA, d_rhs_MRAMLA, d_rhs_ROSEM, d_rhs_RBI, d_rhs_DRAMA, d_rhs_COSEM, d_rhs_ACOSEM,
			d_rhs_MRP_OSEM, d_rhs_Quad_OSEM, d_rhs_L_OSEM, d_rhs_FMH_OSEM, d_rhs_WeightedMean_OSEM, d_rhs_TV_OSEM, d_rhs_AD_OSEM, d_rhs_APLS_OSEM, d_rhs_TGV_OSEM,
			d_rhs_MRP_BSREM, d_rhs_Quad_BSREM, d_rhs_L_BSREM, d_rhs_FMH_BSREM, d_rhs_WeightedMean_BSREM, d_rhs_TV_BSREM, d_rhs_AD_BSREM, d_rhs_APLS_BSREM,
			d_rhs_TGV_BSREM,
			d_rhs_MRP_MBSREM, d_rhs_Quad_MBSREM, d_rhs_L_MBSREM, d_rhs_FMH_MBSREM, d_rhs_WeightedMean_MBSREM, d_rhs_TV_MBSREM, d_rhs_AD_MBSREM, d_rhs_APLS_MBSREM,
			d_rhs_TGV_MBSREM,
			d_rhs_MRP_ROSEM, d_rhs_Quad_ROSEM, d_rhs_L_ROSEM, d_rhs_FMH_ROSEM, d_rhs_WeightedMean_ROSEM, d_rhs_TV_ROSEM, d_rhs_AD_ROSEM, d_rhs_APLS_ROSEM,
			d_rhs_TGV_ROSEM,
			d_rhs_MRP_RBI, d_rhs_Quad_RBI, d_rhs_L_RBI, d_rhs_FMH_RBI, d_rhs_WeightedMean_RBI, d_rhs_TV_RBI, d_rhs_AD_RBI, d_rhs_APLS_RBI, d_rhs_TGV_RBI,
			d_rhs_MRP_COSEM, d_rhs_Quad_COSEM, d_rhs_L_COSEM, d_rhs_FMH_COSEM, d_rhs_WeightedMean_COSEM, d_rhs_TV_COSEM, d_rhs_AD_COSEM, d_rhs_APLS_COSEM,
			d_rhs_TGV_COSEM, d_h,
			d_MRP_COSEM, d_Quad_COSEM, d_L_COSEM, d_FMH_COSEM, d_WeightedMean_COSEM, d_TV_COSEM, d_AD_COSEM, d_APLS_COSEM, d_TGV_COSEM, d_COSEM, d_ACOSEM);
		atomicAdd_g_f(&Summ[local_ind], local_ele);
	}
}

// Denominator (forward projection), orthogonal distance based ray tracer
inline void orth_distance_denominator_mlem(const int tempi, const uint Nx, const float y_diff, const float x_diff, __constant float* y_center, __constant float* x_center, const float kerroin,
	const float length_, float* temp, const uint temp_ijk, const RecMethodsOpenCL MethodList, 
	const float local_sino, float ax[], const uint d_Ny, const uint d_N, const int tempj,
	const __global float* d_OSEM, const __global float* d_MRP_OSEM, const __global float* d_Quad_OSEM, const __global float* d_L_OSEM, const __global float* d_FMH_OSEM, const __global float* d_WeightedMean_OSEM,
	const __global float* d_TV_OSEM, const __global float* d_AD_OSEM, const __global float* d_APLS_OSEM, const __global float* d_TGV_OSEM) {

	const float diff = kerroin - x_diff * y_center[tempj];
	for (int uu = tempi; uu >= 0; uu--) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		*temp += local_ele;
		if (local_sino > 0.f) {
			denominator_mlem(local_ele, MethodList, ax, local_ind,
				d_OSEM, d_MRP_OSEM, d_Quad_OSEM, d_L_OSEM, d_FMH_OSEM, d_WeightedMean_OSEM, d_TV_OSEM, d_AD_OSEM, d_APLS_OSEM, d_TGV_OSEM);
		}
	}
	for (int uu = tempi + 1; uu < Nx; uu++) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		*temp += local_ele;
		if (local_sino > 0.f) {
			denominator_mlem(local_ele, MethodList, ax, local_ind,
				d_OSEM, d_MRP_OSEM, d_Quad_OSEM, d_L_OSEM, d_FMH_OSEM, d_WeightedMean_OSEM, d_TV_OSEM, d_AD_OSEM, d_APLS_OSEM, d_TGV_OSEM);
		}
	}
}

// RHS, orthogonal distance based ray tracer
inline void orth_distance_rhs_mlem(const int tempi, const uint Nx, const float y_diff, const float x_diff, const float y_center, __constant float* x_center, const float kerroin,
	const float length_, const float temp, const uint temp_ijk, const RecMethodsOpenCL MethodList,
	const float ax[], const uint d_Ny, const uint d_N,
	__global float* d_rhs_OSEM, __global float* d_rhs_MRP_OSEM, __global float* d_rhs_Quad_OSEM, __global float* d_rhs_L_OSEM, __global float* d_rhs_FMH_OSEM, __global float* d_rhs_WeightedMean_OSEM,
	__global float* d_rhs_TV_OSEM, __global float* d_rhs_AD_OSEM, __global float* d_rhs_APLS_OSEM, __global float* d_rhs_TGV_OSEM, __global float* Summ) {

	const float diff = kerroin - x_diff * y_center;
	for (int uu = tempi; uu >= 0; uu--) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		local_ele *= temp;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		atomicAdd_g_f(&Summ[local_ind], local_ele);
		rhs_mlem(MethodList, local_ele, ax, local_ind,
			d_rhs_OSEM, d_rhs_MRP_OSEM, d_rhs_Quad_OSEM, d_rhs_L_OSEM, d_rhs_FMH_OSEM, d_rhs_WeightedMean_OSEM, d_rhs_TV_OSEM, d_rhs_AD_OSEM, d_rhs_APLS_OSEM, d_rhs_TGV_OSEM);
	}
	for (int uu = tempi + 1; uu < Nx; uu++) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		local_ele *= temp;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		rhs_mlem(MethodList, local_ele, ax, local_ind,
			d_rhs_OSEM, d_rhs_MRP_OSEM, d_rhs_Quad_OSEM, d_rhs_L_OSEM, d_rhs_FMH_OSEM, d_rhs_WeightedMean_OSEM, d_rhs_TV_OSEM, d_rhs_AD_OSEM, d_rhs_APLS_OSEM, d_rhs_TGV_OSEM);
		atomicAdd_g_f(&Summ[local_ind], local_ele);
	}
}


// Denominator (forward projection), multi-GPU version, orthogonal distance based ray tracer
inline void orth_distance_multi(const int tempi, const uint Nx, const float y_diff, const float x_diff, const float y_center,
	__constant float* x_center, const float kerroin, const float length_, float* temp, const uint temp_ijk, 
	const float local_sino, float *ax, const uint d_Ny, const uint d_N, const uchar no_norm, const bool RHS, const bool SUMMA, const __global float* d_OSEM, __global float* Summ, __global float* d_rhs_OSEM) {

	const float diff = kerroin - x_diff * y_center;
	for (int uu = tempi; uu >= 0; uu--) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		if (RHS) {
			local_ele *= *temp;
			atomicAdd_g_f(&d_rhs_OSEM[local_ind], (local_ele * *ax));
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
			atomicAdd_g_f(&d_rhs_OSEM[local_ind], (local_ele * *ax));
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
				denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
			}
		}
	}
}

// RHS, multi-GPU version, orthogonal distance based ray tracer
inline void orth_distance_rhs_multi(const int tempi, const uint Nx, const float y_diff, const float x_diff,
	const float y_center, __constant float* x_center, const float kerroin, const float length_, const float temp, const uint temp_ijk,
	const float ax, const uint d_Ny, const uint d_N, const uchar no_norm, __global float* d_rhs_OSEM, __global float* Summ) {

	const float diff = kerroin - x_diff * y_center;
	for (int uu = tempi; uu >= 0; uu--) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		local_ele *= temp;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		//rhs_multi(local_ele, ax, &d_rhs_OSEM[local_ind]);
		atomicAdd_g_f(&d_rhs_OSEM[local_ind], (local_ele * ax));
		if (no_norm == 0u)
			atomicAdd_g_f(&Summ[local_ind], local_ele);
	}
	for (int uu = tempi + 1; uu < convert_int(Nx); uu++) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		local_ele *= temp;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		if (no_norm == 0u)
			atomicAdd_g_f(&Summ[local_ind], local_ele);
		atomicAdd_g_f(&d_rhs_OSEM[local_ind], (local_ele * ax));
		//rhs_multi(local_ele, ax, &d_rhs_OSEM[local_ind]);
	}
}


// Denominator (forward projection), multi-GPU version, orthogonal distance based ray tracer
inline void orth_distance_multi_3D(const int tempi, const uint Nx, const uint Nz, const float y_diff, const float x_diff, const float z_diff, 
	const float y_center, __constant float* x_center, __constant float* z_center, const float crystal_size_z, float* temp, const uint temp_ijk, 
	const float local_sino, float* ax, const uint d_Ny, const uint d_N, const int tempk, const uint Nxy, const float xs, const float ys, const float zs,
	const int dec, const uchar no_norm, const bool RHS, const bool SUMMA, const __global float* d_OSEM, __global float* Summ, __global float* d_rhs_OSEM) {

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
			const uint local_ind = compute_ind_orth_mfree_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
			if (RHS) {
				local_ele *= *temp;
				atomicAdd_g_f(&d_rhs_OSEM[local_ind], (local_ele * *ax));
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
				const uint local_ind = compute_ind_orth_mfree_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
				if (RHS) {
					local_ele *= *temp;
					atomicAdd_g_f(&d_rhs_OSEM[local_ind], (local_ele * *ax));
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
			const uint local_ind = compute_ind_orth_mfree_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
			if (RHS) {
				local_ele *= *temp;
				atomicAdd_g_f(&d_rhs_OSEM[local_ind], (local_ele * *ax));
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
				const uint local_ind = compute_ind_orth_mfree_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
				if (RHS) {
					local_ele *= *temp;
					atomicAdd_g_f(&d_rhs_OSEM[local_ind], (local_ele * *ax));
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

// Orthogonal distance, forward/back projection, orthogonal distance based ray tracer
inline void orth_distance_bpfp(const int tempi, const uint Nx, const float y_diff, const float x_diff, const float y_center, __constant float* x_center,
	const float kerroin, const float length_, const float temp, const uint temp_ijk, const uint d_Ny, const uint d_N,
	__global float* d_Summ, __global float* d_output, const uchar no_norm, float d_rhs_local) {

	const float diff = kerroin - x_diff * y_center;
	for (int uu = tempi; uu >= 0; uu--) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		local_ele *= temp;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		//if (fp) {
		//	*d_output_local += (d_rhs[local_ind] * local_ele);
		//}
		//else
		atomicAdd_g_f(&d_output[local_ind], local_ele * d_rhs_local);
		if (!no_norm)
			atomicAdd_g_f(&d_Summ[local_ind], local_ele);
	}
	for (int uu = tempi + 1; uu < Nx; uu++) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		local_ele *= temp;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		//if (fp) {
		//	*d_output_local += (d_rhs[local_ind] * local_ele);
		//}
		//else
		atomicAdd_g_f(&d_output[local_ind], local_ele * d_rhs_local);
		if (!no_norm)
			atomicAdd_g_f(&d_Summ[local_ind], local_ele);
	}
}


// Denominator, forward/back projection, orthogonal distance based ray tracer
inline void orth_distance_denominator_bpfp(const int tempi, const uint Nx, const float y_diff, const float x_diff, __constant float* y_center,
	__constant float* x_center, const float kerroin, const float length_, float* temp, const uint temp_ijk, 
	const uint d_Ny, const uint d_N, const int tempj, const uchar fp, float* d_output_local, const __global float* d_rhs) {

	const float diff = kerroin - x_diff * y_center[tempi];
	for (int uu = tempi; uu >= 0; uu--) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		*temp += local_ele;
		if (fp) {
			*d_output_local += (d_rhs[local_ind] * local_ele);
		}
	}
	for (int uu = tempi + 1; uu < Nx; uu++) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		*temp += local_ele;
		if (fp) {
			*d_output_local += (d_rhs[local_ind] * local_ele);
		}
	}
}


// Denominator, forward/back projection, orthogonal distance based ray tracer
inline void orth_distance_denominator_bpfp_3D(const int tempi, const uint Nx, const uint Nz, const float y_diff, const float x_diff, const float z_diff,
	const float y_center, __constant float* x_center, __constant float* z_center, const float crystal_size_z, float* temp, const uint temp_ijk, 
	const uint d_Ny, const uint d_N, const int tempk, const uint Nxy, const float xs, const float ys, const float zs, const uchar fp, float* d_output_local, const __global float* d_rhs) {

	int loppu1 = -1;
	int loppu2 = -1;
	for (int uu = tempi; uu >= 0; uu--) {
		for (int zz = tempk; zz >= 0; zz--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
			if (local_ele <= THR) {
				loppu1 = zz;
				break;
			}
			const uint local_ind = compute_ind_orth_mfree_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
			*temp += local_ele;
			if (fp) {
				*d_output_local += (d_rhs[local_ind] * local_ele);
			}
		}
		for (uint zz = convert_uint(tempk) + 1u; zz < Nz; zz++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
			if (local_ele <= THR) {
				loppu2 = zz;
				break;
			}
			const uint local_ind = compute_ind_orth_mfree_3D(convert_uint(uu), temp_ijk, convert_int(zz), d_N, Nxy);
			*temp += local_ele;
			if (fp) {
				*d_output_local += (d_rhs[local_ind] * local_ele);
			}
		}
		if (loppu1 == tempk && loppu2 == (tempk + 1))
			break;
	}
	loppu1 = -1;
	loppu2 = -1;
	for (uint uu = convert_uint(tempi) + 1u; uu < Nx; uu++) {
		for (int zz = tempk; zz >= 0; zz--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
			if (local_ele <= THR) {
				loppu1 = zz;
				break;
			}
			const uint local_ind = compute_ind_orth_mfree_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
			*temp += local_ele;
			if (fp) {
				*d_output_local += (d_rhs[local_ind] * local_ele);
			}
		}
		for (uint zz = convert_uint(tempk) + 1u; zz < Nz; zz++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
			if (local_ele <= THR) {
				loppu2 = zz;
				break;
			}
			const uint local_ind = compute_ind_orth_mfree_3D(convert_uint(uu), temp_ijk, convert_int(zz), d_N, Nxy);
			*temp += local_ele;
			if (fp) {
				*d_output_local += (d_rhs[local_ind] * local_ele);
			}
		}
		if (loppu1 == tempk && loppu2 == (tempk + 1))
			break;
	}
}

// Denominator (forward projection), orthogonal distance based ray tracer
inline void orth_distance_denominator_perpendicular(const float diff1, __constant float* center1, const float kerroin,
	const float length_, float* temp, const RecMethodsOpenCL MethodList, const uint d_attenuation_correction, const uint d_normalization, float ax[],
	const float d_b, const float d, const float d_d1, const uint d_N1, const uint d_N2, const uint z_loop,
	const __global float* d_atten, const float d_norm, const float local_sino, const uint d_N, const uint d_NN,
	const __global float* d_OSEM, const __global float* d_RAMLA, const __global float* d_MRAMLA, const __global float* d_ROSEM, const __global float* d_RBI, const __global float* d_DRAMA,
	const __global float* d_COSEM, const __global float* d_ACOSEM,
	const __global float* d_MRP_OSEM, const __global float* d_Quad_OSEM, const __global float* d_L_OSEM, const __global float* d_FMH_OSEM, const __global float* d_WeightedMean_OSEM,
	const __global float* d_TV_OSEM, const __global float* d_AD_OSEM, const __global float* d_APLS_OSEM, const __global float* d_TGV_OSEM,
	const __global float* d_MRP_BSREM, const __global float* d_Quad_BSREM, const __global float* d_L_BSREM, const __global float* d_FMH_BSREM, const __global float* d_WeightedMean_BSREM,
	const __global float* d_TV_BSREM, const __global float* d_AD_BSREM, const __global float* d_APLS_BSREM, const __global float* d_TGV_BSREM,
	const __global float* d_MRP_MBSREM, const __global float* d_Quad_MBSREM, const __global float* d_L_MBSREM, const __global float* d_FMH_MBSREM, const __global float* d_WeightedMean_MBSREM,
	const __global float* d_TV_MBSREM, const __global float* d_AD_MBSREM, const __global float* d_APLS_MBSREM, const __global float* d_TGV_MBSREM,
	const __global float* d_MRP_ROSEM, const __global float* d_Quad_ROSEM, const __global float* d_L_ROSEM, const __global float* d_FMH_ROSEM, const __global float* d_WeightedMean_ROSEM,
	const __global float* d_TV_ROSEM, const __global float* d_AD_ROSEM, const __global float* d_APLS_ROSEM, const __global float* d_TGV_ROSEM,
	const __global float* d_MRP_RBI, const __global float* d_Quad_RBI, const __global float* d_L_RBI, const __global float* d_FMH_RBI, const __global float* d_WeightedMean_RBI,
	const __global float* d_TV_RBI, const __global float* d_AD_RBI, const __global float* d_APLS_RBI, const __global float* d_TGV_RBI,
	const __global float* d_MRP_COSEM, const __global float* d_Quad_COSEM, const __global float* d_L_COSEM, const __global float* d_FMH_COSEM, const __global float* d_WeightedMean_COSEM,
	const __global float* d_TV_COSEM, const __global float* d_AD_COSEM, const __global float* d_APLS_COSEM, const __global float* d_TGV_COSEM) {

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
				denominator(local_ele, MethodList, ax, local_ind,
					d_OSEM, d_RAMLA, d_MRAMLA, d_ROSEM, d_RBI, d_DRAMA, d_COSEM, d_ACOSEM,
					d_MRP_OSEM, d_Quad_OSEM, d_L_OSEM, d_FMH_OSEM, d_WeightedMean_OSEM,	d_TV_OSEM, d_AD_OSEM, d_APLS_OSEM, d_TGV_OSEM,
					d_MRP_BSREM, d_Quad_BSREM, d_L_BSREM, d_FMH_BSREM, d_WeightedMean_BSREM, d_TV_BSREM, d_AD_BSREM, d_APLS_BSREM, d_TGV_BSREM,
					d_MRP_MBSREM, d_Quad_MBSREM, d_L_MBSREM, d_FMH_MBSREM, d_WeightedMean_MBSREM, d_TV_MBSREM, d_AD_MBSREM, d_APLS_MBSREM, d_TGV_MBSREM,
					d_MRP_ROSEM, d_Quad_ROSEM, d_L_ROSEM, d_FMH_ROSEM, d_WeightedMean_ROSEM, d_TV_ROSEM, d_AD_ROSEM, d_APLS_ROSEM, d_TGV_ROSEM,
					d_MRP_RBI, d_Quad_RBI, d_L_RBI, d_FMH_RBI, d_WeightedMean_RBI, d_TV_RBI, d_AD_RBI, d_APLS_RBI, d_TGV_RBI,
					d_MRP_COSEM, d_Quad_COSEM, d_L_COSEM, d_FMH_COSEM, d_WeightedMean_COSEM, d_TV_COSEM, d_AD_COSEM, d_APLS_COSEM, d_TGV_COSEM);
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
				denominator(local_ele, MethodList, ax, local_ind,
					d_OSEM, d_RAMLA, d_MRAMLA, d_ROSEM, d_RBI, d_DRAMA, d_COSEM, d_ACOSEM,
					d_MRP_OSEM, d_Quad_OSEM, d_L_OSEM, d_FMH_OSEM, d_WeightedMean_OSEM, d_TV_OSEM, d_AD_OSEM, d_APLS_OSEM, d_TGV_OSEM,
					d_MRP_BSREM, d_Quad_BSREM, d_L_BSREM, d_FMH_BSREM, d_WeightedMean_BSREM, d_TV_BSREM, d_AD_BSREM, d_APLS_BSREM, d_TGV_BSREM,
					d_MRP_MBSREM, d_Quad_MBSREM, d_L_MBSREM, d_FMH_MBSREM, d_WeightedMean_MBSREM, d_TV_MBSREM, d_AD_MBSREM, d_APLS_MBSREM, d_TGV_MBSREM,
					d_MRP_ROSEM, d_Quad_ROSEM, d_L_ROSEM, d_FMH_ROSEM, d_WeightedMean_ROSEM, d_TV_ROSEM, d_AD_ROSEM, d_APLS_ROSEM, d_TGV_ROSEM,
					d_MRP_RBI, d_Quad_RBI, d_L_RBI, d_FMH_RBI, d_WeightedMean_RBI, d_TV_RBI, d_AD_RBI, d_APLS_RBI, d_TGV_RBI,
					d_MRP_COSEM, d_Quad_COSEM, d_L_COSEM, d_FMH_COSEM, d_WeightedMean_COSEM, d_TV_COSEM, d_AD_COSEM, d_APLS_COSEM, d_TGV_COSEM);
			}
			local_ind += d_NN;
		}
	}
	*temp = 1.f / *temp;
	if (d_attenuation_correction)
		* temp *= jelppi;
	if (d_normalization == 1u)
		* temp *= d_norm;
}


// RHS, orthogonal distance based ray tracer
inline void orth_distance_rhs_perpendicular(const float diff1, __constant float* center1, const float kerroin,
	const float length_, const float temp, const RecMethodsOpenCL MethodList, float ax[],
	const float d_b, const float d, const float d_d1, const uint d_N1, const uint d_N2, const uint z_loop, const uint d_N, const uint d_NN,
	__global float* d_rhs_OSEM, __global float* d_rhs_RAMLA, __global float* d_rhs_MRAMLA, __global float* d_rhs_ROSEM, __global float* d_rhs_RBI, __global float* d_rhs_DRAMA,
	__global float* d_rhs_COSEM, __global float* d_rhs_ACOSEM,
	__global float* d_rhs_MRP_OSEM, __global float* d_rhs_Quad_OSEM, __global float* d_rhs_L_OSEM, __global float* d_rhs_FMH_OSEM, __global float* d_rhs_WeightedMean_OSEM,
	__global float* d_rhs_TV_OSEM, __global float* d_rhs_AD_OSEM, __global float* d_rhs_APLS_OSEM, __global float* d_rhs_TGV_OSEM,
	__global float* d_rhs_MRP_BSREM, __global float* d_rhs_Quad_BSREM, __global float* d_rhs_L_BSREM, __global float* d_rhs_FMH_BSREM, __global float* d_rhs_WeightedMean_BSREM,
	__global float* d_rhs_TV_BSREM, __global float* d_rhs_AD_BSREM, __global float* d_rhs_APLS_BSREM, __global float* d_rhs_TGV_BSREM,
	__global float* d_rhs_MRP_MBSREM, __global float* d_rhs_Quad_MBSREM, __global float* d_rhs_L_MBSREM, __global float* d_rhs_FMH_MBSREM, __global float* d_rhs_WeightedMean_MBSREM,
	__global float* d_rhs_TV_MBSREM, __global float* d_rhs_AD_MBSREM, __global float* d_rhs_APLS_MBSREM, __global float* d_rhs_TGV_MBSREM,
	__global float* d_rhs_MRP_ROSEM, __global float* d_rhs_Quad_ROSEM, __global float* d_rhs_L_ROSEM, __global float* d_rhs_FMH_ROSEM, __global float* d_rhs_WeightedMean_ROSEM,
	__global float* d_rhs_TV_ROSEM, __global float* d_rhs_AD_ROSEM, __global float* d_rhs_APLS_ROSEM, __global float* d_rhs_TGV_ROSEM,
	__global float* d_rhs_MRP_RBI, __global float* d_rhs_Quad_RBI, __global float* d_rhs_L_RBI, __global float* d_rhs_FMH_RBI, __global float* d_rhs_WeightedMean_RBI,
	__global float* d_rhs_TV_RBI, __global float* d_rhs_AD_RBI, __global float* d_rhs_APLS_RBI, __global float* d_rhs_TGV_RBI,
	__global float* d_rhs_MRP_COSEM, __global float* d_rhs_Quad_COSEM, __global float* d_rhs_L_COSEM, __global float* d_rhs_FMH_COSEM, __global float* d_rhs_WeightedMean_COSEM,
	__global float* d_rhs_TV_COSEM, __global float* d_rhs_AD_COSEM, __global float* d_rhs_APLS_COSEM, __global float* d_rhs_TGV_COSEM,
	const float d_h,
	const __global float* d_MRP_COSEM, const __global float* d_Quad_COSEM, const __global float* d_L_COSEM, const __global float* d_FMH_COSEM, const __global float* d_WeightedMean_COSEM,
	const __global float* d_TV_COSEM, const __global float* d_AD_COSEM, const __global float* d_APLS_COSEM, const __global float* d_TGV_COSEM, const __global float* d_COSEM,
	const __global float* d_ACOSEM,
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
			//atomicAdd_g_f(&d_rhs_OSEM[local_ind], 1.f);
			atomicAdd_g_f(&Summ[local_ind], local_ele);
			rhs(MethodList, local_ele, ax, local_ind,
				d_rhs_OSEM, d_rhs_RAMLA, d_rhs_MRAMLA, d_rhs_ROSEM, d_rhs_RBI, d_rhs_DRAMA, d_rhs_COSEM, d_rhs_ACOSEM,
				d_rhs_MRP_OSEM, d_rhs_Quad_OSEM, d_rhs_L_OSEM, d_rhs_FMH_OSEM, d_rhs_WeightedMean_OSEM, d_rhs_TV_OSEM, d_rhs_AD_OSEM, d_rhs_APLS_OSEM, d_rhs_TGV_OSEM,
				d_rhs_MRP_BSREM, d_rhs_Quad_BSREM, d_rhs_L_BSREM, d_rhs_FMH_BSREM, d_rhs_WeightedMean_BSREM, d_rhs_TV_BSREM, d_rhs_AD_BSREM, d_rhs_APLS_BSREM,
				d_rhs_TGV_BSREM,
				d_rhs_MRP_MBSREM, d_rhs_Quad_MBSREM, d_rhs_L_MBSREM, d_rhs_FMH_MBSREM, d_rhs_WeightedMean_MBSREM, d_rhs_TV_MBSREM, d_rhs_AD_MBSREM, d_rhs_APLS_MBSREM,
				d_rhs_TGV_MBSREM,
				d_rhs_MRP_ROSEM, d_rhs_Quad_ROSEM, d_rhs_L_ROSEM, d_rhs_FMH_ROSEM, d_rhs_WeightedMean_ROSEM, d_rhs_TV_ROSEM, d_rhs_AD_ROSEM, d_rhs_APLS_ROSEM,
				d_rhs_TGV_ROSEM,
				d_rhs_MRP_RBI, d_rhs_Quad_RBI, d_rhs_L_RBI, d_rhs_FMH_RBI, d_rhs_WeightedMean_RBI, d_rhs_TV_RBI, d_rhs_AD_RBI, d_rhs_APLS_RBI, d_rhs_TGV_RBI,
				d_rhs_MRP_COSEM, d_rhs_Quad_COSEM, d_rhs_L_COSEM, d_rhs_FMH_COSEM, d_rhs_WeightedMean_COSEM, d_rhs_TV_COSEM, d_rhs_AD_COSEM, d_rhs_APLS_COSEM,
				d_rhs_TGV_COSEM, d_h,
				d_MRP_COSEM, d_Quad_COSEM, d_L_COSEM, d_FMH_COSEM, d_WeightedMean_COSEM, d_TV_COSEM, d_AD_COSEM, d_APLS_COSEM, d_TGV_COSEM, d_COSEM, d_ACOSEM);
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
			rhs(MethodList, local_ele, ax, local_ind,
				d_rhs_OSEM, d_rhs_RAMLA, d_rhs_MRAMLA, d_rhs_ROSEM, d_rhs_RBI, d_rhs_DRAMA, d_rhs_COSEM, d_rhs_ACOSEM,
				d_rhs_MRP_OSEM, d_rhs_Quad_OSEM, d_rhs_L_OSEM, d_rhs_FMH_OSEM, d_rhs_WeightedMean_OSEM, d_rhs_TV_OSEM, d_rhs_AD_OSEM, d_rhs_APLS_OSEM, d_rhs_TGV_OSEM,
				d_rhs_MRP_BSREM, d_rhs_Quad_BSREM, d_rhs_L_BSREM, d_rhs_FMH_BSREM, d_rhs_WeightedMean_BSREM, d_rhs_TV_BSREM, d_rhs_AD_BSREM, d_rhs_APLS_BSREM,
				d_rhs_TGV_BSREM,
				d_rhs_MRP_MBSREM, d_rhs_Quad_MBSREM, d_rhs_L_MBSREM, d_rhs_FMH_MBSREM, d_rhs_WeightedMean_MBSREM, d_rhs_TV_MBSREM, d_rhs_AD_MBSREM, d_rhs_APLS_MBSREM,
				d_rhs_TGV_MBSREM,
				d_rhs_MRP_ROSEM, d_rhs_Quad_ROSEM, d_rhs_L_ROSEM, d_rhs_FMH_ROSEM, d_rhs_WeightedMean_ROSEM, d_rhs_TV_ROSEM, d_rhs_AD_ROSEM, d_rhs_APLS_ROSEM,
				d_rhs_TGV_ROSEM,
				d_rhs_MRP_RBI, d_rhs_Quad_RBI, d_rhs_L_RBI, d_rhs_FMH_RBI, d_rhs_WeightedMean_RBI, d_rhs_TV_RBI, d_rhs_AD_RBI, d_rhs_APLS_RBI, d_rhs_TGV_RBI,
				d_rhs_MRP_COSEM, d_rhs_Quad_COSEM, d_rhs_L_COSEM, d_rhs_FMH_COSEM, d_rhs_WeightedMean_COSEM, d_rhs_TV_COSEM, d_rhs_AD_COSEM, d_rhs_APLS_COSEM,
				d_rhs_TGV_COSEM, d_h,
				d_MRP_COSEM, d_Quad_COSEM, d_L_COSEM, d_FMH_COSEM, d_WeightedMean_COSEM, d_TV_COSEM, d_AD_COSEM, d_APLS_COSEM, d_TGV_COSEM, d_COSEM, d_ACOSEM);
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

// Denominator (forward projection), orthogonal distance based ray tracer
inline void orth_distance_denominator_perpendicular_mlem(const float diff1, __constant float* center1, const float kerroin,
	const float length_, float* temp, const RecMethodsOpenCL MethodList, const uint d_attenuation_correction, const uint d_normalization, float ax[],
	const float d_b, const float d, const float d_d1, const uint d_N1, const uint d_N2, const uint z_loop,
	const __global float* d_atten, const float d_norm, const float local_sino, const uint d_N, const uint d_NN,
	const __global float* d_OSEM, const __global float* d_MRP_OSEM, const __global float* d_Quad_OSEM, const __global float* d_L_OSEM, 
	const __global float* d_FMH_OSEM, const __global float* d_WeightedMean_OSEM,
	const __global float* d_TV_OSEM, const __global float* d_AD_OSEM, const __global float* d_APLS_OSEM, const __global float* d_TGV_OSEM) {

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
				denominator_mlem(local_ele, MethodList, ax, local_ind,
					d_OSEM, d_MRP_OSEM, d_Quad_OSEM, d_L_OSEM, d_FMH_OSEM, d_WeightedMean_OSEM, d_TV_OSEM, d_AD_OSEM, d_APLS_OSEM, d_TGV_OSEM);
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
		if (local_sino > 0.f) {
			for (uint kk = 0u; kk < d_N2; kk++) {
				denominator_mlem(local_ele, MethodList, ax, local_ind,
					d_OSEM, d_MRP_OSEM, d_Quad_OSEM, d_L_OSEM, d_FMH_OSEM, d_WeightedMean_OSEM, d_TV_OSEM, d_AD_OSEM, d_APLS_OSEM, d_TGV_OSEM);
				local_ind += d_NN;
			}
		}
	}
	*temp = 1.f / *temp;
	if (d_attenuation_correction)
		* temp *= jelppi;
	if (d_normalization == 1u)
		* temp *= d_norm;
}


// RHS, orthogonal distance based ray tracer
inline void orth_distance_rhs_perpendicular_mlem(const float diff1, __constant float* center1, const float kerroin,
	const float length_, const float temp, const RecMethodsOpenCL MethodList, float ax[],
	const float d_b, const float d, const float d_d1, const uint d_N1, const uint d_N2, const uint z_loop, const uint d_N, const uint d_NN,
	__global float* d_rhs_OSEM, __global float* d_rhs_MRP_OSEM, __global float* d_rhs_Quad_OSEM, __global float* d_rhs_L_OSEM, __global float* d_rhs_FMH_OSEM, __global float* d_rhs_WeightedMean_OSEM,
	__global float* d_rhs_TV_OSEM, __global float* d_rhs_AD_OSEM, __global float* d_rhs_APLS_OSEM, __global float* d_rhs_TGV_OSEM, __global float* Summ) {

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
			atomicAdd_g_f(&Summ[local_ind], local_ele);
			rhs_mlem(MethodList, local_ele, ax, local_ind,
				d_rhs_OSEM, d_rhs_MRP_OSEM, d_rhs_Quad_OSEM, d_rhs_L_OSEM, d_rhs_FMH_OSEM, d_rhs_WeightedMean_OSEM, d_rhs_TV_OSEM, d_rhs_AD_OSEM, d_rhs_APLS_OSEM, d_rhs_TGV_OSEM);
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
			rhs_mlem(MethodList, local_ele, ax, local_ind,
				d_rhs_OSEM, d_rhs_MRP_OSEM, d_rhs_Quad_OSEM, d_rhs_L_OSEM, d_rhs_FMH_OSEM, d_rhs_WeightedMean_OSEM, d_rhs_TV_OSEM, d_rhs_AD_OSEM, d_rhs_APLS_OSEM, d_rhs_TGV_OSEM);
			atomicAdd_g_f(&Summ[local_ind], local_ele);
			local_ind += d_NN;
		}
	}
}


// Denominator (forward projection), orthogonal distance based ray tracer, multi-GPU
inline void orth_distance_denominator_perpendicular_multi(const float diff1, __constant float* center1, const float kerroin,
	const float length_, float* temp, const uint d_attenuation_correction, const uint d_normalization, float* ax, const float d_b, const float d, const float d_d1,
	const uint d_N1, const uint d_N2, const uint z_loop, const __global float* d_atten, const float d_norm, const float local_sino, const uint d_N, const uint d_NN,
	const __global float* d_OSEM) {

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
				denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
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
		if (local_sino > 0.f) {
			for (uint kk = 0u; kk < d_N2; kk++) {
				denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
				local_ind += d_NN;
			}
		}
	}
	*temp = 1.f / *temp;
	if (d_attenuation_correction)
		*temp *= jelppi;
	if (d_normalization == 1u)
		*temp *= d_norm;
}


// RHS, orthogonal distance based ray tracer, multi-GPU
inline void orth_distance_rhs_perpendicular_multi(const float diff1, __constant float* center1, const float kerroin,
	const float length_, const float temp, const float ax, const float d_b, const float d, const float d_d1, const uint d_N1, const uint d_N2,
	const uint z_loop, const uint d_N, const uint d_NN, const uchar no_norm, __global float* d_rhs_OSEM, __global float* Summ) {

	const uint zz = z_loop * d_N2 * d_N1;
	const int apu = perpendicular_start(d_b, d, d_d1, d_N1);
	for (int uu = apu; uu >= 0; uu--) {
		float local_ele = compute_element_orth(kerroin, diff1, center1[uu], length_);
		if (local_ele <= THR)
			break;
		local_ele *= temp;
		uint local_ind = convert_uint_sat(uu) * d_N + zz;
		for (uint kk = 0u; kk < d_N2; kk++) {
			if (no_norm == 0u)
				atomicAdd_g_f(&Summ[local_ind], local_ele);
			atomicAdd_g_f(&d_rhs_OSEM[local_ind], (local_ele * ax));
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
			atomicAdd_g_f(&d_rhs_OSEM[local_ind], (local_ele * ax));
			if (no_norm == 0u)
				atomicAdd_g_f(&Summ[local_ind], local_ele);
			local_ind += d_NN;
		}
	}
}

// Denominator (forward projection), orthogonal distance based ray tracer, multi-GPU
inline void orth_distance_denominator_perpendicular_bpfp(const float diff1, __constant float* center1, const float kerroin,
	const float length_, float* temp, const uint d_attenuation_correction, const uint d_normalization, const float d_b, const float d, const float d_d1,
	const uint d_N1, const uint d_N2, const uint z_loop, const __global float* d_atten, const float d_norm, const uint d_N, const uint d_NN,
	const __global float* d_rhs, const uchar fp, float* d_output_local) {

	const uint zz = z_loop * d_N2 * d_N1;
	const int apu = perpendicular_start(d_b, d, d_d1, d_N1);
	float jelppi = 0.f;
	for (int uu = apu; uu >= 0; uu--) {
		float local_ele = compute_element_orth(kerroin, diff1, center1[uu], length_);
		if (local_ele <= THR)
			break;
		*temp += (local_ele * d_N2);
		if ((d_attenuation_correction && uu == apu) || fp) {
			uint local_ind = convert_uint_sat(uu) * d_N + zz;
			for (uint kk = 0u; kk < d_N2; kk++) {
				if (d_attenuation_correction && uu == apu)
					jelppi += (d_d1 * -d_atten[local_ind]);
				if (fp)
					*d_output_local += (d_rhs[local_ind] * local_ele);
				local_ind += d_NN;
			}
		}
	}
	for (uint uu = convert_uint_sat(apu) + 1u; uu < d_N1; uu++) {
		float local_ele = compute_element_orth(kerroin, diff1, center1[uu], length_);
		if (local_ele <= THR)
			break;
		*temp += (local_ele * d_N2);
		if (fp) {
			uint local_ind = uu * d_N + zz;
			for (uint kk = 0u; kk < d_N2; kk++) {
				*d_output_local += (d_rhs[local_ind] * local_ele);
				local_ind += d_NN;
			}
		}
	}
	*temp = 1.f / *temp;
	if (d_attenuation_correction)
		* temp *= jelppi;
	if (d_normalization == 1u)
		* temp *= d_norm;
}

// Denominator (forward projection), orthogonal distance based ray tracer, multi-GPU
inline void orth_distance_denominator_perpendicular_bpfp_3D(__constant float* center1, const float center2, __constant float* z_center,
	float* temp, const uint d_attenuation_correction, const uint d_normalization, const float d_b, const float d, const float d_d1,
	const uint d_N1, const uint d_N2, const uint z_loop, const __global float* d_atten, const float d_norm, const float local_sino, const uint d_N, const uint d_NN, const __global float* d_rhs, const uchar fp, float* d_output_local,
	const float xs, const float ys, const float zs, const float xl, const float yl, const float zl, const float crystal_size_z, const uint Nyx, const uint Nz) {

	const uint apu = perpendicular_start(d_b, d, d_d1, d_N1);
	float jelppi = 0.f;
	for (int zz = convert_int_sat(z_loop); zz >= 0; zz--) {
		for (int uu = convert_int_sat(apu); uu >= 0; uu--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			*temp += (local_ele * d_N2);
			if ((d_attenuation_correction && zz == convert_int_sat(z_loop) && uu == convert_int_sat(apu)) || fp) {
				uint local_ind = uu * d_N + zz * Nyx;
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (d_attenuation_correction && zz == convert_int_sat(z_loop) && uu == convert_int_sat(apu))
						jelppi += (d_d1 * -d_atten[local_ind]);
					if (fp)
						* d_output_local += (d_rhs[local_ind] * local_ele);
					local_ind += d_NN;
				}
			}
		}
		for (uint uu = convert_uint_sat(apu) + 1u; uu < d_N1; uu++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			*temp += (local_ele * d_N2);
			if (fp) {
				uint local_ind = uu * d_N + zz * Nyx;
				for (uint kk = 0u; kk < d_N2; kk++) {
					* d_output_local += (d_rhs[local_ind] * local_ele);
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
			*temp += (local_ele * d_N2);
			if (fp) {
				uint local_ind = uu * d_N + zz * Nyx;
				for (uint kk = 0u; kk < d_N2; kk++) {
					* d_output_local += (d_rhs[local_ind] * local_ele);
					local_ind += d_NN;
				}
			}
		}
		for (uint uu = convert_uint_sat(apu) + 1u; uu < d_N1; uu++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			*temp += (local_ele * d_N2);
			if (fp) {
				uint local_ind = uu * d_N + zz * Nyx;
				for (uint kk = 0u; kk < d_N2; kk++) {
					* d_output_local += (d_rhs[local_ind] * local_ele);
					local_ind += d_NN;
				}
			}
		}
	}
	*temp = 1.f / *temp;
	if (d_attenuation_correction)
		* temp *= jelppi;
	if (d_normalization == 1u)
		* temp *= d_norm;
}

inline void orth_distance_denominator_perpendicular_multi_3D(__constant float* center1, const float center2, __constant float* z_center,
	float* temp, const uint d_attenuation_correction, const uint d_normalization, float* ax, const float d_b, const float d, const float d_d1,
	const uint d_N1, const uint d_N2, const uint z_loop, const __global float* d_atten, const float d_norm, const float local_sino, const uint d_N, const uint d_NN, const __global float* d_OSEM, 
	const float xs, const float ys, const float zs, const float xl, const float yl, const float zl, const float crystal_size_z, const uint Nyx, const uint Nz) {

	//const uint zz = z_loop * d_N2 * d_N1;
	const uint apu = perpendicular_start(d_b, d, d_d1, d_N1);
	float jelppi = 0.f;
	for (int zz = convert_int_sat(z_loop); zz >= 0; zz--) {
		for (int uu = convert_int_sat(apu); uu >= 0; uu--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			*temp += (local_ele * d_N2);
			uint local_ind = uu * d_N + zz * Nyx;
			for (uint kk = 0u; kk < d_N2; kk++) {
				if (d_attenuation_correction && zz == convert_int_sat(z_loop) && uu == convert_int_sat(apu))
					jelppi += (d_d1 * -d_atten[local_ind]);
				if (local_sino > 0.f) {
					denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
				}
				local_ind += d_NN;
			}
		}
		for (uint uu = apu + 1; uu < d_N1; uu++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			*temp += (local_ele * d_N2);
			uint local_ind = uu * d_N + zz * Nyx;
			if (local_sino > 0.f) {
				for (uint kk = 0u; kk < d_N2; kk++) {
					denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
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
			*temp += (local_ele * d_N2);
			uint local_ind = uu * d_N + zz * Nyx;
			if (local_sino > 0.f) {
				for (uint kk = 0u; kk < d_N2; kk++) {
					denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
					local_ind += d_NN;
				}
			}
		}
		for (uint uu = apu + 1; uu < d_N1; uu++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			*temp += (local_ele * d_N2);
			uint local_ind = uu * d_N + zz * Nyx;
			if (local_sino > 0.f) {
				for (uint kk = 0u; kk < d_N2; kk++) {
					denominator_multi(local_ele, ax, &d_OSEM[local_ind]);
					local_ind += d_NN;
				}
			}
		}
	}
	*temp = 1. / *temp;
	if (d_attenuation_correction)
		*temp *= jelppi;
	if (d_normalization == 1u)
		*temp *= d_norm;
}


inline void orth_distance_rhs_perpendicular_multi_3D(__constant float* center1, const float center2, __constant float* z_center,
	const float temp, const float ax, const float d_b, const float d, const float d_d1, const uint d_N1, const uint d_N2,
	const uint z_loop, const uint d_N, const uint d_NN, const float xs, const float ys, const float zs, const float xl, const float yl,
	const float zl, const float crystal_size_z, const uint Nyx, const uint Nz, const uchar no_norm, __global float* Summ, __global float* d_rhs_OSEM) {

	//const uint zz = z_loop * d_N2 * d_N1;
	const int apu = perpendicular_start(d_b, d, d_d1, d_N1);
	for (int zz = convert_int(z_loop); zz >= 0; zz--) {
		for (int uu = apu; uu >= 0; uu--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			local_ele *= temp;
			for (uint kk = 0u; kk < d_N2; kk++) {
				if (no_norm == 0u)
					atomicAdd_g_f(&Summ[local_ind], local_ele);
				atomicAdd_g_f(&d_rhs_OSEM[local_ind], (local_ele * ax));
				local_ind += d_NN;
			}
		}
		for (uint uu = convert_uint(apu) + 1u; uu < d_N1; uu++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			local_ele *= temp;
			for (uint kk = 0u; kk < d_N2; kk++) {
				atomicAdd_g_f(&d_rhs_OSEM[local_ind], (local_ele * ax));
				if (no_norm == 0u)
					atomicAdd_g_f(&Summ[local_ind], local_ele);
				local_ind += d_NN;
			}
		}
	}
	for (uint zz = z_loop + 1u; zz < Nz; zz++) {
		for (int uu = apu; uu >= 0; uu--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			local_ele *= temp;
			for (uint kk = 0u; kk < d_N2; kk++) {
				atomicAdd_g_f(&d_rhs_OSEM[local_ind], (local_ele * ax));
				if (no_norm == 0u)
					atomicAdd_g_f(&Summ[local_ind], local_ele);
				local_ind += d_NN;
			}
		}
		for (uint uu = convert_uint(apu) + 1u; uu < d_N1; uu++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			local_ele *= temp;
			for (uint kk = 0u; kk < d_N2; kk++) {
				if (no_norm == 0u)
					atomicAdd_g_f(&Summ[local_ind], local_ele);
				atomicAdd_g_f(&d_rhs_OSEM[local_ind], (local_ele * ax));
				local_ind += d_NN;
			}
		}
	}
}

void orth_distance_summ_perpendicular_3D(__constant float* center1, const float center2, __constant float* z_center, const float temp,
	const float d_b, const float d, const float d_d1, const uint d_N1, const uint d_N2, const uint z_loop, const uint d_N, const uint d_NN,
	const float xs, const float ys, const float zs, const float xl, const float yl, const float zl, const float crystal_size_z, const uint Nyx, const uint Nz, __global float* Summ) {

	//const uint zz = z_loop * d_N2 * d_N1;
	const uint apu = perpendicular_start(d_b, d, d_d1, d_N1);
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

// RHS, orthogonal distance based ray tracer, multi-GPU
inline void orth_distance_rhs_perpendicular_bpfp(const float diff1, __constant float* center1, const float kerroin,
	const float length_, const float temp, float* ax, const float d_b, const float d, const float d_d1, const uint d_N1, const uint d_N2,
	const uint z_loop, const uint d_N, const uint d_NN, __global float* d_output, __global float* Summ, const uchar no_norm, const float d_rhs_local) {

	const uint zz = z_loop * d_N2 * d_N1;
	const int apu = perpendicular_start(d_b, d, d_d1, d_N1);
	for (int uu = apu; uu >= 0; uu--) {
		float local_ele = compute_element_orth(kerroin, diff1, center1[uu], length_);
		if (local_ele <= THR)
			break;
		local_ele *= temp;
		uint local_ind = convert_uint_sat(uu) * d_N + zz;
		for (uint kk = 0u; kk < d_N2; kk++) {
			atomicAdd_g_f(&d_output[local_ind], local_ele * d_rhs_local);
			if (!no_norm)
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
			atomicAdd_g_f(&d_output[local_ind], local_ele * d_rhs_local);
			if (!no_norm)
				atomicAdd_g_f(&Summ[local_ind], local_ele);
			local_ind += d_NN;
		}
	}
}

// RHS, orthogonal distance based ray tracer, multi-GPU
inline void orth_distance_rhs_perpendicular_bpfp_3D(__constant float* center1, const float center2, __constant float* z_center,
	const float temp, const float d_b, const float d, const float d_d1, const uint d_N1, const uint d_N2,
	const uint z_loop, const uint d_N, const uint d_NN, const float xs, const float ys, const float zs, const float xl, const float yl,
	const float zl, const float crystal_size_z, const uint Nyx, const uint Nz, __global float* d_output, __global float* Summ, const uchar no_norm, const float d_rhs_local) {

	const int apu = perpendicular_start(d_b, d, d_d1, d_N1);
	for (int zz = convert_int(z_loop); zz >= 0; zz--) {
		for (int uu = apu; uu >= 0; uu--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			local_ele *= temp;
			for (uint kk = 0u; kk < d_N2; kk++) {
				atomicAdd_g_f(&d_output[local_ind], local_ele * d_rhs_local);
				if (!no_norm)
					atomicAdd_g_f(&Summ[local_ind], local_ele);
				local_ind += d_NN;
			}
		}
		for (uint uu = convert_uint_sat(apu) + 1u; uu < d_N1; uu++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			local_ele *= temp;
			for (uint kk = 0u; kk < d_N2; kk++) {
				atomicAdd_g_f(&d_output[local_ind], local_ele * d_rhs_local);
				if (!no_norm)
					atomicAdd_g_f(&Summ[local_ind], local_ele);
				local_ind += d_NN;
			}
		}
	}
	for (uint zz = z_loop + 1u; zz < Nz; zz++) {
		for (int uu = apu; uu >= 0; uu--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			local_ele *= temp;
			for (uint kk = 0u; kk < d_N2; kk++) {
				atomicAdd_g_f(&d_output[local_ind], local_ele * d_rhs_local);
				if (!no_norm)
					atomicAdd_g_f(&Summ[local_ind], local_ele);
				local_ind += d_NN;
			}
		}
		for (uint uu = convert_uint_sat(apu) + 1u; uu < d_N1; uu++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz * Nyx;
			local_ele *= temp;
			for (uint kk = 0u; kk < d_N2; kk++) {
				atomicAdd_g_f(&d_output[local_ind], local_ele * d_rhs_local);
				if (!no_norm)
					atomicAdd_g_f(&Summ[local_ind], local_ele);
				local_ind += d_NN;
			}
		}
	}
}

inline void orth_element_per1(const float diff1, __constant float* center1, const float kerroin, const float length_, float* temp, const uint d_attenuation_correction, const uint d_normalization, 
	const float d_b, const float d, const float d_d1, const uint d_N1, const uint d_N2, const uint z_loop, const __global float* d_atten, const float d_norm, const uint d_N, 
	const uint d_NN, float* axCOSEM, const __global float* d_COSEM, const RecMethodsOpenCL MethodListOpenCL, const int apu) {

	float jelppi = 0.f;
	const uint zz = z_loop * d_N1 * d_N2;
	for (int uu = apu; uu >= 0; uu--) {
		float local_ele = compute_element_orth(kerroin, diff1, center1[uu], length_);
		if (local_ele <= THR)
			break;
		*temp += (local_ele * d_N2);
		uint local_ind = convert_uint_sat(uu) * d_N + zz;
		if ((d_attenuation_correction && uu == apu) || MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) {
			for (uint kk = 0u; kk < d_N2; kk++) {
				if (d_attenuation_correction && uu == apu)
					jelppi += (d_d1 * -d_atten[local_ind]);
				if (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0)
					*axCOSEM += (local_ele * d_COSEM[local_ind]);
				local_ind += d_NN;
			}
		}
	}
	for (uint uu = convert_uint_sat(apu) + 1u; uu < d_N1; uu++) {
		float local_ele = compute_element_orth(kerroin, diff1, center1[uu], length_);
		if (local_ele <= THR)
			break;
		*temp += (local_ele * d_N2);
		uint local_ind = uu * d_N + zz;
		if (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) {
			for (uint kk = 0u; kk < d_N2; kk++) {
				if (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0)
					*axCOSEM += (local_ele * d_COSEM[local_ind]);
				local_ind += d_NN;
			}
		}
	}
	*temp = 1.f / *temp;
	if (d_attenuation_correction)
		*temp *= jelppi;
	if (d_normalization == 1u)
		* temp *= d_norm;
}

inline void orth_element_per_3D(__constant float* center1, const float center2, __constant float* z_center, float* temp, const uint d_attenuation_correction, const uint d_normalization,
	const float d_b, const float d, const float d_d1, const uint d_N1, const uint d_N2, const uint z_loop, const __global float* d_atten, const float d_norm, const uint d_N,
	const uint d_NN, float* axCOSEM, const __global float* d_COSEM, const RecMethodsOpenCL MethodListOpenCL, const int apu, const float xs, const float ys, const float zs, const float xl, const float yl, const float zl, const float crystal_size_z, const uint Nyx, const uint Nz) {

	float jelppi = 0.f;
	for (int zz = convert_int_sat(z_loop); zz >= 0; zz--) {
		for (int uu = apu; uu >= 0; uu--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			*temp += (local_ele * d_N2);
			uint local_ind = convert_uint_sat(uu) * d_N + zz;
			if (d_attenuation_correction || MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) {
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (d_attenuation_correction &&  zz == convert_int_sat(z_loop) && uu == convert_int_sat(apu))
						jelppi += (d_d1 * -d_atten[local_ind]);
					if (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0)
						* axCOSEM += (local_ele * d_COSEM[local_ind]);
					local_ind += d_NN;
				}
			}
		}
		for (uint uu = convert_uint_sat(apu) + 1u; uu < d_N1; uu++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			*temp += (local_ele * d_N2);
			uint local_ind = uu * d_N + zz;
			if ( MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) {
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0)
						* axCOSEM += (local_ele * d_COSEM[local_ind]);
					local_ind += d_NN;
				}
			}
		}
	}
	for (uint zz = z_loop + 1u; zz < Nz; zz++) {
		for (int uu = apu; uu >= 0; uu--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			*temp += (local_ele * d_N2);
			uint local_ind = convert_uint_sat(uu) * d_N + zz;
			if (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) {
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0)
						* axCOSEM += (local_ele * d_COSEM[local_ind]);
					local_ind += d_NN;
				}
			}
		}
		for (uint uu = convert_uint_sat(apu) + 1u; uu < d_N1; uu++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			*temp += (local_ele * d_N2);
			uint local_ind = uu * d_N + zz;
			if (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) {
				for (uint kk = 0u; kk < d_N2; kk++) {
					if (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0)
						* axCOSEM += (local_ele * d_COSEM[local_ind]);
					local_ind += d_NN;
				}
			}
		}
	}
	*temp = 1.f / *temp;
	if (d_attenuation_correction)
		* temp *= jelppi;
	if (d_normalization == 1u)
		* temp *= d_norm;
}

inline void orth_distance_perpendicular_cosem_3D(__constant float* center1, const float center2, __constant float* z_center, const uint d_attenuation_correction, const uint d_normalization, 
	const float d_b, const float d, const float d_d1, const uint d_N1, const uint d_N2, const uint z_loop, const __global float* d_atten, const float d_norm, const uint d_N, const uint Np, 
	const float d_epps, const float local_sino, const uint idx, const uint d_NN, const float xs, const float ys, const float zs, const float xl, const float yl, const float zl, const float crystal_size_z, const uint Nyx, const uint Nz, const RecMethodsOpenCL MethodListOpenCL, const uint d_alku,
	const uchar MBSREM_prepass, const __global float* d_COSEM, const float d_h, __global float* d_aco, __global float* d_co, __global float* d_Summ, __global float* d_Amin, 
	__global float* d_ACOSEM_lhs) {

	float temp = 0.f;
	float axCOSEM = 0.f;
	float axACOSEM = 0.f;
	const int apu = perpendicular_start(d_b, d, d_d1, d_N1);
	orth_element_per_3D(center1, center2, z_center, &temp, d_attenuation_correction, d_normalization, d_b, d, d_d1, d_N1, d_N2, z_loop, d_atten, d_norm, d_N, d_NN, &axCOSEM, d_COSEM, MethodListOpenCL, apu, xs, ys, zs, xl, yl, zl, crystal_size_z, Nyx, Nz);
	float minimi = 1e8f;
	if (d_alku == 0 && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && local_sino != 0.f) {
		if (axCOSEM == 0.f)
			axCOSEM = d_epps;
		else
			axCOSEM *= temp;
		axCOSEM = local_sino / axCOSEM;
	}
	for (int zz = convert_int_sat(z_loop); zz >= 0; zz--) {
		for (int uu = apu; uu >= 0; uu--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = convert_uint_sat(uu) * d_N + zz;
			local_ele *= temp;
			if (d_alku == 0 && (MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1) {
				if (local_ele < minimi && local_ele != 0.f)
					minimi = local_ele;
			}
			for (uint ii = 0u; ii < Np; ii++) {
				if (d_alku == 0) {
					if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino != 0.f)
						atomicAdd_g_f(&d_co[local_ind], (axCOSEM * local_ele * d_COSEM[local_ind]));
					if (MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1)
						atomicAdd_g_f(&d_aco[local_ind], axCOSEM * (native_powr(d_COSEM[local_ind], d_h) * local_ele));
					if (MBSREM_prepass == 1)
						atomicAdd_g_f(&d_Summ[local_ind], local_ele);
				}
				else
					axACOSEM += (local_ele * d_COSEM[local_ind++]);
				local_ind += d_NN;
			}
		}
		for (uint uu = convert_uint_sat(apu) + 1u; uu < d_N1; uu++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz;
			local_ele *= temp;
			if ((MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1) {
				if (local_ele < minimi && local_ele != 0.f)
					minimi = local_ele;
			}
			for (uint ii = 0u; ii < Np; ii++) {
				if (d_alku == 0) {
					if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino != 0.f)
						atomicAdd_g_f(&d_co[local_ind], (axCOSEM * local_ele * d_COSEM[local_ind]));
					if (MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1)
						atomicAdd_g_f(&d_aco[local_ind], axCOSEM * (native_powr(d_COSEM[local_ind], d_h) * local_ele));
					if (MBSREM_prepass == 1)
						atomicAdd_g_f(&d_Summ[local_ind], local_ele);
				}
				else
					axACOSEM += (local_ele * d_COSEM[local_ind++]);
				local_ind += d_NN;
			}
		}
	}
	for (uint zz = z_loop + 1u; zz < Nz; zz++) {
		for (int uu = apu; uu >= 0; uu--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = convert_uint_sat(uu) * d_N + zz;
			local_ele *= temp;
			if (d_alku == 0 && (MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1) {
				if (local_ele < minimi && local_ele != 0.f)
					minimi = local_ele;
			}
			for (uint ii = 0u; ii < Np; ii++) {
				if (d_alku == 0) {
					if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino != 0.f)
						atomicAdd_g_f(&d_co[local_ind], (axCOSEM * local_ele * d_COSEM[local_ind]));
					if (MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1)
						atomicAdd_g_f(&d_aco[local_ind], axCOSEM * (native_powr(d_COSEM[local_ind], d_h) * local_ele));
					if (MBSREM_prepass == 1)
						atomicAdd_g_f(&d_Summ[local_ind], local_ele);
				}
				else
					axACOSEM += (local_ele * d_COSEM[local_ind++]);
				local_ind += d_NN;
			}
		}
		for (uint uu = convert_uint_sat(apu) + 1u; uu < d_N1; uu++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint local_ind = uu * d_N + zz;
			local_ele *= temp;
			if ((MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1) {
				if (local_ele < minimi && local_ele != 0.f)
					minimi = local_ele;
			}
			for (uint ii = 0u; ii < Np; ii++) {
				if (d_alku == 0) {
					if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino != 0.f)
						atomicAdd_g_f(&d_co[local_ind], (axCOSEM * local_ele * d_COSEM[local_ind]));
					if (MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1)
						atomicAdd_g_f(&d_aco[local_ind], axCOSEM * (native_powr(d_COSEM[local_ind], d_h) * local_ele));
					if (MBSREM_prepass == 1)
						atomicAdd_g_f(&d_Summ[local_ind], local_ele);
				}
				else
					axACOSEM += (local_ele * d_COSEM[local_ind++]);
				local_ind += d_NN;
			}
		}
	}
	if ((MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1)
		d_Amin[idx] = minimi;
	if (d_alku == 0)
		d_ACOSEM_lhs[idx] = axACOSEM;
}

inline void orth_distance_perpendicular_cosem(const float diff1, __constant float* center1, const float kerroin, const float length_, const uint d_attenuation_correction, const uint d_normalization,
	const float d_b, const float d, const float d_d1, const uint d_N1, const uint d_N2, const uint z_loop, const __global float* d_atten, const float d_norm, const uint d_N, const uint Np,
	const float d_epps, const float local_sino, const uint idx, const uint d_NN, const RecMethodsOpenCL MethodListOpenCL, const uint d_alku,
	const uchar MBSREM_prepass, const __global float* d_COSEM, const float d_h, __global float* d_aco, __global float* d_co, __global float* d_Summ, __global float* d_Amin,
	__global float* d_ACOSEM_lhs) {

	float temp = 0.f;
	float axCOSEM = 0.f;
	float axACOSEM = 0.f;
	const uint zz = z_loop * d_N2 * d_N1;
	const int apu = perpendicular_start(d_b, d, d_d1, d_N1);
	orth_element_per1(diff1, center1, kerroin, length_, &temp, d_attenuation_correction, d_normalization, d_b, d, d_d1, d_N1, d_N2, z_loop, d_atten, d_norm, d_N, d_NN, &axCOSEM, d_COSEM, MethodListOpenCL, apu);
	float minimi = 1e8f;
	if (d_alku == 0 && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && local_sino != 0.f) {
		if (axCOSEM == 0.f)
			axCOSEM = d_epps;
		else
			axCOSEM *= temp;
		axCOSEM = local_sino / axCOSEM;
	}
	for (int uu = apu; uu >= 0; uu--) {
		float local_ele = compute_element_orth(kerroin, diff1, center1[uu], length_);
		if (local_ele <= THR)
			break;
		uint local_ind = convert_uint_sat(uu) * d_N + zz;
		local_ele *= temp;
		if (d_alku == 0 && (MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1) {
			if (local_ele < minimi && local_ele != 0.f)
				minimi = local_ele;
		}
		for (uint ii = 0u; ii < Np; ii++) {
			if (d_alku == 0) {
				if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino != 0.f)
					atomicAdd_g_f(&d_co[local_ind], (axCOSEM * local_ele * d_COSEM[local_ind]));
				if (MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1)
					atomicAdd_g_f(&d_aco[local_ind], axCOSEM * (native_powr(d_COSEM[local_ind], d_h) * local_ele));
				if (MBSREM_prepass == 1)
					atomicAdd_g_f(&d_Summ[local_ind], local_ele);
			}
			else
				axACOSEM += (local_ele * d_COSEM[local_ind++]);
			local_ind += d_NN;
		}
	}
	for (uint uu = convert_uint_sat(apu) + 1u; uu < d_N1; uu++) {
		float local_ele = compute_element_orth(kerroin, diff1, center1[uu], length_);
		if (local_ele <= THR)
			break;
		uint local_ind = uu * d_N + zz;
		local_ele *= temp;
		if ((MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1) {
			if (local_ele < minimi && local_ele != 0.f)
				minimi = local_ele;
		}
		for (uint ii = 0u; ii < Np; ii++) {
			if (d_alku == 0) {
				if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino != 0.f)
					atomicAdd_g_f(&d_co[local_ind], (axCOSEM * local_ele * d_COSEM[local_ind]));
				if (MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1)
					atomicAdd_g_f(&d_aco[local_ind], axCOSEM * (native_powr(d_COSEM[local_ind], d_h) * local_ele));
				if (MBSREM_prepass == 1)
					atomicAdd_g_f(&d_Summ[local_ind], local_ele);
			}
			else
				axACOSEM += (local_ele * d_COSEM[local_ind++]);
			local_ind += d_NN;
		}
	}
	if ((MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1)
		d_Amin[idx] = minimi;
	if (d_alku == 0)
		d_ACOSEM_lhs[idx] = axACOSEM;
}

// Denominator (forward projection), orthogonal distance based ray tracer
inline void orth_distance_denominator_cosem(const int tempi, const uint Nx, const float y_diff, const float x_diff, __constant float* y_center, __constant float* x_center, const float kerroin,
	const float length_, float* temp, const uint temp_ijk, const RecMethodsOpenCL MethodListOpenCL, const uint d_alku,
	const float local_sino, const uint d_Ny, const uint d_N, const int tempj, float* axCOSEM, 
	const __global float* d_COSEM) {

	const float diff = kerroin - x_diff * y_center[tempj];
	for (int uu = tempi; uu >= 0; uu--) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		*temp += local_ele;
		if (local_sino > 0.f && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0) {
			*axCOSEM += (local_ele * d_COSEM[local_ind]);
		}
	}
	for (int uu = tempi + 1; uu < Nx; uu++) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		*temp += local_ele;
		if (local_sino > 0.f) {
			*axCOSEM += (local_ele * d_COSEM[local_ind]);
		}
	}
}

inline void orth_distance_rhs_cosem(const int tempi, const uint Nx, const float y_diff, const float x_diff, const float y_center, __constant float* x_center,
	const float kerroin, const float length_, const float temp, const uint temp_ijk, float* axACOSEM, const uint d_Ny, const uint d_N, const RecMethodsOpenCL MethodListOpenCL,
	__global float* d_Summ, float* minimi, const uchar MBSREM_prepass, const __global float* d_COSEM, __global float* d_co, __global float* d_aco, 
	const uint d_alku, const float d_epps, const float local_sino, const float d_h, const float axCOSEM) {

	const float diff = kerroin - x_diff * y_center;
	for (int uu = tempi; uu >= 0; uu--) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		local_ele *= temp;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		if (d_alku == 0 && (MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1) {
			if (local_ele < *minimi && local_ele != 0.f)
				*minimi = local_ele;
		}
		if (d_alku == 0) {
			if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
				atomicAdd_g_f(&d_co[local_ind], (axCOSEM * local_ele * d_COSEM[local_ind]));
			if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
				atomicAdd_g_f(&d_aco[local_ind], axCOSEM * (native_powr(d_COSEM[local_ind], d_h) * local_ele));
			if (MBSREM_prepass == 1)
				atomicAdd_g_f(&d_Summ[local_ind], local_ele);
		}
		else
			*axACOSEM += (local_ele * d_COSEM[local_ind]);
	}
	for (int uu = tempi + 1; uu < convert_int(Nx); uu++) {
		float local_ele = compute_element_orth(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		local_ele *= temp;
		const uint local_ind = compute_ind_orth(uu, temp_ijk, d_Ny, Nx, d_N);
		if ((MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1) {
			if (local_ele < *minimi && local_ele != 0.f)
				*minimi = local_ele;
		}
		if (d_alku == 0) {
			if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
				atomicAdd_g_f(&d_co[local_ind], (axCOSEM * local_ele * d_COSEM[local_ind]));
			if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
				atomicAdd_g_f(&d_aco[local_ind], axCOSEM * (native_powr(d_COSEM[local_ind], d_h) * local_ele));
			if (MBSREM_prepass == 1)
				atomicAdd_g_f(&d_Summ[local_ind], local_ele);
		}
		else
			*axACOSEM += (local_ele * d_COSEM[local_ind]);
	}
}

// Denominator (forward projection), orthogonal distance based ray tracer
inline void orth_distance_denominator_cosem_3D(const int tempi, const uint Nx, const uint Nz, const float y_diff, const float x_diff, const float z_diff,
	const float y_center, __constant float* x_center, __constant float* z_center, const float crystal_size_z, float* temp, const uint temp_ijk, const RecMethodsOpenCL MethodListOpenCL, const uint d_alku,
	const float local_sino, const uint d_Ny, const uint d_N, const int tempk, const uint Nxy, const float xs, const float ys, const float zs, float* axCOSEM,
	const __global float* d_COSEM) {

	int loppu1 = -1;
	int loppu2 = -1;
	for (int uu = tempi; uu >= 0; uu--) {
		for (int zz = tempk; zz >= 0; zz--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
			if (local_ele <= THR) {
				loppu1 = zz;
				break;
			}
			const uint local_ind = compute_ind_orth_mfree_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
			*temp += local_ele;
			if (local_sino > 0.f && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0) {
				*axCOSEM += (local_ele * d_COSEM[local_ind]);
			}
		}
		for (uint zz = convert_uint(tempk) + 1u; zz < Nz; zz++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
			if (local_ele <= THR) {
				loppu1 = zz;
				break;
			}
			const uint local_ind = compute_ind_orth_mfree_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
			*temp += local_ele;
			if (local_sino > 0.f && (MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0) {
				*axCOSEM += (local_ele * d_COSEM[local_ind]);
			}
		}
	}
	loppu1 = -1;
	loppu2 = -1;
	for (uint uu = convert_uint(tempi) + 1u; uu < Nx; uu++) {
		for (int zz = tempk; zz >= 0; zz--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
			if (local_ele <= THR) {
				loppu1 = zz;
				break;
			}
			const uint local_ind = compute_ind_orth_mfree_3D(uu, temp_ijk, convert_uint(zz), d_N, Nxy);
			*temp += local_ele;
			if (local_sino > 0.f) {
				*axCOSEM += (local_ele * d_COSEM[local_ind]);
			}
		}
		for (uint zz = convert_uint(tempk) + 1u; zz < Nz; zz++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
			if (local_ele <= THR) {
				loppu1 = zz;
				break;
			}
			const uint local_ind = compute_ind_orth_mfree_3D(uu, temp_ijk, convert_uint(zz), d_N, Nxy);
			*temp += local_ele;
			if (local_sino > 0.f) {
				*axCOSEM += (local_ele * d_COSEM[local_ind]);
			}
		}
	}
}

inline void orth_distance_rhs_cosem_3D(const int tempi, const uint Nx, const uint Nz, const float y_diff, const float x_diff, const float z_diff,
	const float y_center, __constant float* x_center, __constant float* z_center, const float crystal_size_z, const float temp, const uint temp_ijk, float* axACOSEM, const uint d_Ny, const uint d_N, const uint Nxy, const float xs, const float ys, const float zs, const int tempk, const RecMethodsOpenCL MethodListOpenCL,
	__global float* d_Summ, float* minimi, const uchar MBSREM_prepass, const __global float* d_COSEM, __global float* d_co, __global float* d_aco,
	const uint d_alku, const float d_epps, const float local_sino, const float d_h, const float axCOSEM) {

	int loppu1 = -1;
	int loppu2 = -1;
	for (int uu = tempi; uu >= 0; uu--) {
		for (int zz = tempk; zz >= 0; zz--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
			if (local_ele <= THR) {
				loppu1 = zz;
				break;
			}
			const uint local_ind = compute_ind_orth_mfree_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
			local_ele *= temp;
			if (d_alku == 0 && (MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1) {
				if (local_ele < *minimi && local_ele != 0.f)
					* minimi = local_ele;
			}
			if (d_alku == 0) {
				if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
					atomicAdd_g_f(&d_co[local_ind], (axCOSEM * local_ele * d_COSEM[local_ind]));
				if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
					atomicAdd_g_f(&d_aco[local_ind], axCOSEM * (native_powr(d_COSEM[local_ind], d_h) * local_ele));
				if (MBSREM_prepass == 1)
					atomicAdd_g_f(&d_Summ[local_ind], local_ele);
			}
			else
				*axACOSEM += (local_ele * d_COSEM[local_ind]);
		}
		for (uint zz = convert_uint(tempk) + 1u; zz < Nz; zz++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
			if (local_ele <= THR) {
				loppu1 = zz;
				break;
			}
			const uint local_ind = compute_ind_orth_mfree_3D(convert_uint(uu), temp_ijk, convert_uint(zz), d_N, Nxy);
			local_ele *= temp;
			if (d_alku == 0 && (MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1) {
				if (local_ele < *minimi && local_ele != 0.f)
					* minimi = local_ele;
			}
			if (d_alku == 0) {
				if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
					atomicAdd_g_f(&d_co[local_ind], (axCOSEM * local_ele * d_COSEM[local_ind]));
				if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
					atomicAdd_g_f(&d_aco[local_ind], axCOSEM * (native_powr(d_COSEM[local_ind], d_h) * local_ele));
				if (MBSREM_prepass == 1)
					atomicAdd_g_f(&d_Summ[local_ind], local_ele);
			}
			else
				*axACOSEM += (local_ele * d_COSEM[local_ind]);
		}
	}
	loppu1 = -1;
	loppu2 = -1;
	for (uint uu = convert_uint(tempi) + 1u; uu < Nx; uu++) {
		for (int zz = tempk; zz >= 0; zz--) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
			if (local_ele <= THR) {
				loppu1 = zz;
				break;
			}
			const uint local_ind = compute_ind_orth_mfree_3D(uu, temp_ijk, convert_uint(zz), d_N, Nxy);
			local_ele *= temp;
			if ((MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1) {
				if (local_ele < *minimi && local_ele != 0.f)
					* minimi = local_ele;
			}
			if (d_alku == 0) {
				if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
					atomicAdd_g_f(&d_co[local_ind], (axCOSEM * local_ele * d_COSEM[local_ind]));
				if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
					atomicAdd_g_f(&d_aco[local_ind], axCOSEM * (native_powr(d_COSEM[local_ind], d_h) * local_ele));
				if (MBSREM_prepass == 1)
					atomicAdd_g_f(&d_Summ[local_ind], local_ele);
			}
			else
				*axACOSEM += (local_ele * d_COSEM[local_ind]);
		}
		for (uint zz = convert_uint(tempk) + 1u; zz < Nz; zz++) {
			float local_ele = compute_element_orth_3D(xs, ys, zs, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
			if (local_ele <= THR) {
				loppu1 = zz;
				break;
			}
			const uint local_ind = compute_ind_orth_mfree_3D(uu, temp_ijk, convert_uint(zz), d_N, Nxy);
			local_ele *= temp;
			if ((MethodListOpenCL.MRAMLA == 1 || MethodListOpenCL.MBSREM == 1) && MBSREM_prepass == 1) {
				if (local_ele < *minimi && local_ele != 0.f)
					* minimi = local_ele;
			}
			if (d_alku == 0) {
				if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
					atomicAdd_g_f(&d_co[local_ind], (axCOSEM * local_ele * d_COSEM[local_ind]));
				if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
					atomicAdd_g_f(&d_aco[local_ind], axCOSEM * (native_powr(d_COSEM[local_ind], d_h) * local_ele));
				if (MBSREM_prepass == 1)
					atomicAdd_g_f(&d_Summ[local_ind], local_ele);
			}
			else
				*axACOSEM += (local_ele * d_COSEM[local_ind]);
		}
	}
}

inline bool siddon_pre_loop_2D(const float b1, const float b2, const float diff1, const float diff2, const float max1, const float max2,
	const float d1, const float d2, const uint N1, const uint N2, int* temp1, int* temp2, float* t1u, float* t2u, uint* Np,
	const int TYPE, const float ys, const float xs, const float yd, const float xd, float* tc, int* u1, int* u2, float* t10, float* t20) {
	// If neither x- nor y-directions are perpendicular
// Correspond to the equations (9) and (10) from reference [2]
	const float apu_tx = b1 - xs;
	const float apu_ty = b2 - ys;
	*t10 = (apu_tx) / (diff1);
	*t20 = (apu_ty) / (diff2);
	const float txback = (max1 - xs) / (diff1);
	const float tyback = (max2 - ys) / (diff2);

	// Equations (5-8)
	const float txmin = fmin(*t10, txback);
	const float txmax = fmax(*t10, txback);
	const float tymin = fmin(*t20, tyback);
	const float tymax = fmax(*t20, tyback);

	// (3-4)
	*tc = fmax(txmin, tymin);
	const float tmax = fmin(txmax, tymax);

	uint imin, imax, jmin, jmax;

	if (TYPE == 0) {
		// If true, then the ray/LOR does not intersect the pixel space --> continue to the next LOR
		if (*tc >= tmax) {
			return true;
		}

		// (11-14)
		if (xs < xd)
			d_g_s_precomp(*tc, txmin, tmax, txmax, &imin, &imax, t10, u1, diff1, b1, d1, xs, N1);
		// (15-18)
		else
			s_g_d_precomp(*tc, txmin, tmax, txmax, &imin, &imax, t10, u1, diff1, b1, d1, xs, N1);

		//Same as above
		if (ys < yd)
			d_g_s_precomp(*tc, tymin, tmax, tymax, &jmin, &jmax, t20, u2, diff2, b2, d2, ys, N2);
		else
			s_g_d_precomp(*tc, tymin, tmax, tymax, &jmin, &jmax, t20, u2, diff2, b2, d2, ys, N2);

		*Np = imax + 1u + jmax + 1u - imin - jmin;
	}
	else {
		// (11-14)
		if (xs < xd)
			d_g_s(*tc, txmin, &imin, t10, u1, diff1, b1, d1, xs);
		// (15-18)
		else
			s_g_d(*tc, txmin, &imax, t10, u1, diff1, b1, d1, xs, N1);
		//Same as above
		if (ys < yd)
			d_g_s(*tc, tymin, &jmin, t20, u2, diff2, b2, d2, ys);
		else
			s_g_d(*tc, tymin, &jmax, t20, u2, diff2, b2, d2, ys, N2);
	}

	// (2) and (19)
	const float pt = ((fmin(*t10, *t20) + *tc) / 2.f);

	// (26)
	*temp1 = voxel_index(pt, diff1, d1, apu_tx);
	// (27)
	*temp2 = voxel_index(pt, diff2, d2, apu_ty);

	// (28)
	*t1u = d1 / fabs(diff1);
	*t2u = d2 / fabs(diff2);

	if (TYPE == 0) {
		if (*temp1 < 0 || *temp2 < 0 || *temp1 >= N1 || *temp2 >= N2)
			return true;
	}

	if (*tc == *t10 || *tc == *t20)
		*tc -= 1e-7f;

	return false;
}

bool siddon_pre_loop_3D(const float bx, const float by, const float bz, const float x_diff, const float y_diff, const float z_diff,
	const float maxxx, const float maxyy, const float bzb, const float dx, const float dy, const float dz,
	const uint Nx, const uint Ny, const uint Nz, int* tempi, int* tempj, int* tempk, float* tyu, float* txu, float* tzu,
	uint* Np, const int TYPE, const float ys, const float xs, const float yd, const float xd, const float zs, const float zd, float* tc, 
	int* iu, int* ju, int* ku, float* tx0, float* ty0, float* tz0) {

	const float apu_tx = bx - xs;
	const float apu_ty = by - ys;
	const float apu_tz = bz - zs;
	*tx0 = (apu_tx) / (x_diff);
	*ty0 = (apu_ty) / (y_diff);
	*tz0 = (apu_tz) / (z_diff);
	const float txback = (maxxx - xs) / (x_diff);
	const float tyback = (maxyy - ys) / (y_diff);
	const float tzback = (bzb - zs) / (z_diff);

	const float txmin = fmin(*tx0, txback);
	const float txmax = fmax(*tx0, txback);
	const float tymin = fmin(*ty0, tyback);
	const float tymax = fmax(*ty0, tyback);
	const float tzmin = fmin(*tz0, tzback);
	const float tzmax = fmax(*tz0, tzback);

	*tc = fmax(fmax(txmin, tzmin), tymin);
	const float tmax = fmin(fmin(txmax, tzmax), tymax);

	uint imin, imax, jmin, jmax, kmin, kmax;

	if (TYPE == 0) {
		if (*tc >= tmax) {
			return true;
		}

		if (xs < xd)
			d_g_s_precomp(*tc, txmin, tmax, txmax, &imin, &imax, tx0, iu, x_diff, bx, dx, xs, Nx);
		else
			s_g_d_precomp(*tc, txmin, tmax, txmax, &imin, &imax, tx0, iu, x_diff, bx, dx, xs, Nx);

		if (ys < yd)
			d_g_s_precomp(*tc, tymin, tmax, tymax, &jmin, &jmax, ty0, ju, y_diff, by, dy, ys, Ny);
		else
			s_g_d_precomp(*tc, tymin, tmax, tymax, &jmin, &jmax, ty0, ju, y_diff, by, dy, ys, Ny);

		if (zs < zd)
			d_g_s_precomp(*tc, tzmin, tmax, tzmax, &kmin, &kmax, tz0, ku, z_diff, bz, dz, zs, Nz);
		else
			s_g_d_precomp(*tc, tzmin, tmax, tzmax, &kmin, &kmax, tz0, ku, z_diff, bz, dz, zs, Nz);

		*Np = (kmax - kmin + 1) + (jmax - jmin + 1) + (imax - imin + 1);
	}
	else {
		if (xs < xd)
			d_g_s(*tc, txmin, &imin, tx0, iu, x_diff, bx, dx, xs);
		else
			s_g_d(*tc, txmin, &imax, tx0, iu, x_diff, bx, dx, xs, Nx);

		if (ys < yd)
			d_g_s(*tc, tymin, &jmin, ty0, ju, y_diff, by, dy, ys);
		else
			s_g_d(*tc, tymin, &jmax, ty0, ju, y_diff, by, dy, ys, Ny);

		if (zs < zd)
			d_g_s(*tc, tzmin, &kmin, tz0, ku, z_diff, bz, dz, zs);
		else
			s_g_d(*tc, tzmin, &kmax, tz0, ku, z_diff, bz, dz, zs, Nz);
	}

	const float pt = ((fmin(fmin(*tz0, *ty0), *tx0) + *tc) / 2.f);

	*tempi = voxel_index(pt, x_diff, dx, apu_tx);
	*tempj = voxel_index(pt, y_diff, dy, apu_ty);
	*tempk = voxel_index(pt, z_diff, dz, apu_tz);

	if (TYPE == 0) {
		if (*tempi < 0 || *tempj < 0 || *tempk < 0 || *tempi >= Nx || *tempj >= Ny || *tempk >= Nz)
			return true;
	}

	*txu = dx / fabs(x_diff);
	*tyu = dy / fabs(y_diff);
	*tzu = dz / fabs(z_diff);

	return false;
}
