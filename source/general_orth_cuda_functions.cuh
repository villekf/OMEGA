/**************************************************************************
* General functions for orthogonal/volume-based ray-tracer CUDA kernel
* files.
*
* Copyright(C) 2020 Ville-Veikko Wettenhovi
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

// Compute the Euclidean norm of a vector
__device__ float e_norm(const float x, const float y, const float z) {
	return sqrtf(x * x + y * y + z * z);
}

// Compute the linear weight for the current voxel
__device__ float compute_element_orth_3D_per(const float xs, const float ys, const float zs, const float xl, const float yl, const float zl, const float crystal_size_z,
	const float xp, const float yp, const float zp) {

	//float x1, y1, z1, x0, y0, z0;

	const float x0 = xp - xs;
	const float y0 = yp - ys;
	const float z0 = zp - zs;

	// Cross product
	const float x1 = yl * z0 - zl * y0;
	const float y1 = zl * x0 - xl * z0;
	const float z1 = xl * y0 - yl * x0;

	const float normi = e_norm(x1, y1, z1);
#ifndef VOL
	return (1.f - normi / crystal_size_z);
#else
	return (normi / crystal_size_z);
#endif
	//return x0;
}
__device__ float compute_element_orth_3D(const float xs, const float ys, const float zs, const float xl, const float yl, const float zl, const float crystal_size_z,
	const float xp) {

	//float x1, y1, z1, x0, y0, z0;

	const float x0 = xp - xs;
	//const float y0 = yp - ys;
	//const float z0 = zp - zs;

	// Cross product
	//const float x1 = yl * z0 - zl * y0;
	const float y1 = zl * x0 - xl;
	const float z1 = ys - yl * x0;

	const float normi = e_norm(zs, y1, z1);
#ifndef VOL
	return (1.f - normi / crystal_size_z);
#else
	return (normi / crystal_size_z);
#endif
	//return x0;
}

// Gaussian weight
//__device__ float compute_element_orth_3D(const float xs, const float ys, const float zs, const float xl, const float yl, const float zl, const float crystal_size_z,
//	const float xp, const float yp, const float zp) {
//
//	float x1, y1, z1, x0, y0, z0;
//
//	x0 = xp - xs;
//	y0 = yp - ys;
//	z0 = zp - zs;
//
//	// Cross product
//	x1 = yl * z0 - zl * y0;
//	y1 = zl * x0 - xl * z0;
//	z1 = xl * y0 - yl * x0;
//
//	const float normi = e_norm(x1, y1, z1);
//
//	float gauss = 0.f;
//
//	if (normi < crystal_size_z)
//		gauss = (1.f - expf(-(normi * normi) / (2.f * (crystal_size_z / 2.35482f) * (crystal_size_z / 2.35482f))));
//
//	return gauss;
//}

__device__ unsigned int compute_ind_orth_3D(const unsigned int tempi, const unsigned int tempijk, const int tempk, const unsigned int d_N, const unsigned int Nyx) {
	unsigned int local_ind = tempi * d_N + tempijk + (unsigned int)(tempk) * Nyx;
	return local_ind;
}

// compute orthogonal distance (2D)
//__device__ float compute_element_orth(const float x_diff, const float y_diff, const float x_center, const float length_) {
//	float local_ele = 1.f - fabs(x_diff + y_diff * x_center) / length_;
//	return local_ele;
//}

// compute voxel index, orthogonal distance based ray tracer
__device__ unsigned int compute_ind_orth(const int tempi, const unsigned int temp_ijk, const unsigned int d_N) {
	unsigned int local_ind = (unsigned int)(tempi) * d_N + temp_ijk;
	return local_ind;
}

#ifndef MBSREM
__device__ void computeIndicesOrth_af(const bool RHS, const bool SUMMA, float local_ele, float* temp, float* ax, const bool no_norm,
	CAST* Summ, CAST* d_rhs_OSEM, const float local_sino, const float* d_OSEM, const unsigned int local_ind,
	const unsigned int im_dim, const unsigned char* MethodList) {
	if (RHS) {
		local_ele *= *temp;
		rhs(MethodList, local_ele, ax, local_ind, im_dim, d_rhs_OSEM);
		if (no_norm == 0u)
#ifdef ATOMIC
			atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
			atomicAdd(&Summ[local_ind], local_ele);
#endif
	}
	else if (SUMMA) {
		local_ele *= *temp;
#ifdef ATOMIC
		atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
		atomicAdd(&Summ[local_ind], local_ele);
#endif
	}
	else {
		*temp += local_ele;
		if (local_sino > 0.f) {
			denominator(local_ele, ax, local_ind, im_dim, d_OSEM);
		}
	}
}
#else
__device__ void computeIndicesOrth_cosem(const bool RHS, float local_ele, float* temp, float* axACOSEM, CAST* Summ, const float local_sino,
	const float* d_COSEM, const unsigned int local_ind, float* d_E, const RecMethodsOpenCL MethodListOpenCL, 
	CAST* d_co, CAST* d_aco, float* minimi, const unsigned int d_alku, const unsigned char MBSREM_prepass, float* axCOSEM, const unsigned int idx) {
	if (RHS) {
		local_ele *= *temp;
		if (d_alku == 0 && (MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
			if (local_ele < *minimi && local_ele > 0.f)
				*minimi = local_ele;
			d_E[idx] += local_ele;
		}
		if (d_alku == 0) {
			if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino > 0.f)
#ifdef ATOMIC
				atomicAdd(&d_co[local_ind], __float2ull_rn(*axCOSEM * local_ele * TH));
#else
				atomicAdd(&d_co[local_ind], (*axCOSEM * local_ele));
#endif
			if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino > 0.f)
#ifdef ATOMIC
				atomicAdd(&d_aco[local_ind], __float2ull_rn(*axCOSEM * TH * local_ele));
#else
				atomicAdd(&d_aco[local_ind], *axCOSEM * local_ele);
#endif
			if (MBSREM_prepass == 1)
#ifdef ATOMIC
				atomicAdd(&Summ[local_ind], __float2ull_rn(local_ele * TH));
#else
				atomicAdd(&Summ[local_ind], local_ele);
#endif
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
#endif


#ifndef PSF_LIMIT
__device__ void computeOrthVoxelDecreasing(const float s1, const float s2, const float s3, const float diff1, const float diff2, const float diff3,
	const float kerroin, const float* center1, const float center2, const float center3, float* temp, const float local_sino, float* ax,
	bool* breikki, const int ind, const int zz, const int yy, int* uu, const bool no_norm, const bool RHS, const bool SUMMA,
	const unsigned int d_Nxy, const unsigned int d_N3, const unsigned int d_N2, CAST* Summ, const float* d_OSEM,
	const float bmin, const float bmax, const float Vmax, const float* V, float* d_store_elements, unsigned int* d_store_indices,
#ifdef MBSREM
	unsigned int* indeksi, float* d_E, const RecMethodsOpenCL MethodListOpenCL, CAST* d_co, CAST* d_aco, float* minimi, const unsigned int d_alku,
	const unsigned char MBSREM_prepass, float* axCOSEM, const unsigned int idx
#else 
	unsigned int* indeksi, CAST* d_rhs_OSEM, const unsigned int im_dim, const unsigned char* MethodList
#endif
) {
	int xx = 0;
	int incr = 0;
	float prev_local = 1.f;
	for (xx = ind; xx >= 0; xx--) {
		float local_ele = compute_element_orth_3D(s1, s2, s3, diff1, diff2, diff3, kerroin, center1[xx]);
#ifdef VOL
		if (local_ele >= bmax && incr > 0 && prev_local <= local_ele) {
			if (xx == ind - 1) {
				*breikki = true;
			}
			break;
		}
		else if (local_ele >= bmax) {
			incr++;
			prev_local = local_ele;
			continue;
	}
		incr++;
		prev_local = local_ele;
		if (local_ele < bmin)
			local_ele = Vmax;
		else
			local_ele = V[__float2uint_rn((local_ele - bmin) * CC)];
#else
		if (local_ele <= THR && incr > 0 && prev_local > local_ele) {
			if (xx == ind - 1) {
				*breikki = true;
			}
			break;
		}
		else if (local_ele <= THR) {
			incr++;
			prev_local = local_ele;
			continue;
		}
		incr = 1;
		prev_local = local_ele;
#endif
		const unsigned int local_ind = compute_ind_orth_3D((unsigned int)(xx), yy * d_N3, (zz), d_N2, d_Nxy);
#ifdef MBSREM
		computeIndicesOrth_cosem(RHS, local_ele, temp, ax, Summ, local_sino, d_OSEM, local_ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
		computeIndicesOrth_af(RHS, SUMMA, local_ele, temp, ax, no_norm, Summ, d_rhs_OSEM, local_sino, d_OSEM, local_ind, im_dim, MethodList);
#endif
#ifdef DEC
		d_store_elements[*indeksi] = local_ele;
		d_store_indices[*indeksi] = local_ind;
		*indeksi = *indeksi + 1u;
#endif
	}
	*uu = xx;
}

__device__ void computeOrthVoxelIncreasing(const float s1, const float s2, const float s3, const float diff1, const float diff2, const float diff3,
	const float kerroin, const float* center1, const float center2, const float center3, float* temp, const float local_sino, float* ax,
	bool* breikki, const int ind, const int zz, const int yy, int* uu, const bool no_norm, const bool RHS, const bool SUMMA,
	const unsigned int d_Nxy, const unsigned int d_N3, const unsigned int d_N2, const unsigned int NN, CAST* Summ, const float* d_OSEM,
	const float bmin, const float bmax, const float Vmax, const float* V, float* d_store_elements, unsigned int* d_store_indices,
#ifdef MBSREM
	unsigned int* indeksi, float* d_E, const RecMethodsOpenCL MethodListOpenCL, CAST* d_co, CAST* d_aco, float* minimi, const unsigned int d_alku,
	const unsigned char MBSREM_prepass, float* axCOSEM, const unsigned int idx
#else 
	unsigned int* indeksi, CAST* d_rhs_OSEM, const unsigned int im_dim, const unsigned char* MethodList
#endif
) {
	int xx = 0;
	int incr = 0;
	float prev_local = 1.f;
	for (xx = ind; xx < NN; xx++) {
		float local_ele = compute_element_orth_3D(s1, s2, s3, diff1, diff2, diff3, kerroin, center1[xx]);
#ifdef VOL
		if (local_ele >= bmax && incr > 0 && prev_local <= local_ele) {
			if (xx == ind - 1) {
				*breikki = true;
			}
			break;
		}
		else if (local_ele >= bmax) {
			incr++;
			prev_local = local_ele;
			continue;
	}
		incr++;
		prev_local = local_ele;
		if (local_ele < bmin)
			local_ele = Vmax;
		else
			local_ele = V[__float2uint_rn((local_ele - bmin) * CC)];
#else
		if (local_ele <= THR && incr > 0 && prev_local > local_ele) {
			if (xx == ind - 1) {
				*breikki = true;
			}
			break;
		}
		else if (local_ele <= THR) {
			incr++;
			prev_local = local_ele;
			continue;
		}
		incr = 1;
		prev_local = local_ele;
#endif
		const unsigned int local_ind = compute_ind_orth_3D((unsigned int)(xx), yy * d_N3, (zz), d_N2, d_Nxy);
#ifdef MBSREM
		computeIndicesOrth_cosem(RHS, local_ele, temp, ax, Summ, local_sino, d_OSEM, local_ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
		computeIndicesOrth_af(RHS, SUMMA, local_ele, temp, ax, no_norm, Summ, d_rhs_OSEM, local_sino, d_OSEM, local_ind, im_dim, MethodList);
#endif
#ifdef DEC
		d_store_elements[*indeksi] = local_ele;
		d_store_indices[*indeksi] = local_ind;
		*indeksi = *indeksi + 1u;
#endif
	}
	*uu = xx;
}

#else
__device__ void computeOrthVoxelDecreasingPSF(float local_ele, const float* local_psf, float* temp, const int temp1, const int temp2, const int temp3, const float local_sino, float* ax,
	bool* breikki, const int uu, const int yy, const bool no_norm, const bool RHS, const bool SUMMA, const unsigned int im_dim, const unsigned int d_N0, const unsigned int d_N1, const unsigned int d_Nz,
	const unsigned int d_N2, const unsigned int d_N3, const unsigned int d_N4, const int x, const int y, CAST* Summ, const float* d_OSEM,
	float* d_store_elements, unsigned int* d_store_indices,
#ifdef MBSREM
	unsigned int* indeksi, float* d_E, const RecMethodsOpenCL MethodListOpenCL, CAST* d_co, CAST* d_aco, float* minimi, const unsigned int d_alku,
	const unsigned char MBSREM_prepass, float* axCOSEM, const unsigned int idx
#else 
	unsigned int* indeksi, CAST* d_rhs_OSEM, const unsigned char* MethodList
#endif
) {
	//if (uu >= 0 && uu < d_N0 && yy >= 0 && yy < d_N1) {
#pragma unroll
	for (int zz = 0; zz >= -PSF_LIMIT; zz--) {
		const int z = (temp3 + zz);
		if (z >= 0 && z < d_Nz) {
			const int local_ind = uu * (int)(d_N2) + yy * (int)(d_N3) + z * (int)(d_N4);
			local_ele *= local_psf[max(max(abs(zz), abs(x)), abs(y))];
#ifdef MBSREM
			computeIndicesOrth_cosem(RHS, local_ele, temp, ax, Summ, local_sino, d_OSEM, local_ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
			computeIndicesOrth_af(RHS, SUMMA, local_ele, temp, ax, no_norm, Summ, d_rhs_OSEM, local_sino, d_OSEM, local_ind, im_dim, MethodList);
#endif
#ifdef DEC
			//if (*indeksi < DEC) {
			d_store_elements[*indeksi] = local_ele;
			d_store_indices[*indeksi] = local_ind;
			*indeksi = *indeksi + 1u;
			//}
#endif
		}
	}
	//}
}

__device__ void computeOrthVoxelIncreasingPSF(float local_ele, const float* local_psf, float* temp, const int temp1, const int temp2, const int temp3, const float local_sino, float* ax,
	bool* breikki, const int uu, const int yy, const bool no_norm, const bool RHS, const bool SUMMA, const unsigned int im_dim, const unsigned int d_N0, const unsigned int d_N1, const unsigned int d_Nz,
	const unsigned int d_N2, const unsigned int d_N3, const unsigned int d_N4, const unsigned int NN, const int x, const int y, CAST* Summ, const float* d_OSEM,
	float* d_store_elements, unsigned int* d_store_indices,
#ifdef MBSREM
	unsigned int* indeksi, float* d_E, const RecMethodsOpenCL MethodListOpenCL, CAST* d_co, CAST* d_aco, float* minimi, const unsigned int d_alku,
	const unsigned char MBSREM_prepass, float* axCOSEM, const unsigned int idx
#else 
	unsigned int* indeksi, CAST* d_rhs_OSEM, const unsigned char* MethodList
#endif
) {
	//if (uu >= 0 && uu < d_N0 && yy >= 0 && yy < d_N1) {
#pragma unroll
	for (int zz = 0; zz <= PSF_LIMIT; zz++) {
		const int z = (temp3 + zz);
		if (z >= 0 && z < d_Nz) {
			const int local_ind = uu * (int)(d_N2) + yy * (int)(d_N3) + z * (int)(d_N4);
			local_ele *= local_psf[max(max(abs(zz), abs(x)), abs(y))];
#ifdef MBSREM
			computeIndicesOrth_cosem(RHS, local_ele, temp, ax, Summ, local_sino, d_OSEM, local_ind, d_E, MethodListOpenCL, d_co, d_aco, minimi, d_alku, MBSREM_prepass, axCOSEM, idx);
#else
			computeIndicesOrth_af(RHS, SUMMA, local_ele, temp, ax, no_norm, Summ, d_rhs_OSEM, local_sino, d_OSEM, local_ind, im_dim, MethodList);
#endif
#ifdef DEC
			//if (*indeksi < DEC) {
			d_store_elements[*indeksi] = local_ele;
			d_store_indices[*indeksi] = local_ind;
			*indeksi = *indeksi + 1u;
			//}
#endif
		}
	}
	//}
}
#endif