
#ifndef AF

__kernel void summa(const __global CAST* d_Summ_device, __global CAST* d_Summ_local, const __global CAST* d_rhs_device, __global CAST* d_rhs_local,
	const uint im_dim, const uchar no_norm) {

	uint gid = get_global_id(0);

	for (uint i = gid; i < im_dim; i += get_global_size(0)) {
		if (no_norm == 0u)
			d_Summ_local[i] += d_Summ_device[i];
		d_rhs_local[i] += d_rhs_device[i];
	}
}


__kernel void mlem(const __global CAST* d_Summ, const __global CAST* d_rhs, __global float* d_mlem, const uint im_dim, const float d_epps) {

	uint gid = get_global_id(0);

	for (uint i = gid; i < im_dim; i += get_global_size(0)) {
#ifdef LISTMODE2
#if defined(ATOMIC) || defined(ATOMIC32)
		float Summ = convert_float(d_Summ[i]);
		d_mlem[i] = Summ / TH;
#else
		d_mlem[i] = d_Summ[i];
#endif
#else
#if defined(ATOMIC) || defined(ATOMIC32)
		float rhs = convert_float(d_rhs[i]);
		float Summ = convert_float(d_Summ[i]);
		if (rhs != 0.f) {
			if (Summ == 0.f)
				d_mlem[i] = d_mlem[i] / d_epps * (rhs / TH);
			else
				d_mlem[i] = d_mlem[i] / (Summ / TH) * (rhs / TH);
		}
		else {
			if (Summ != 0.f)
				d_mlem[i] = d_mlem[i] / (Summ / TH) * d_epps;
		}
#else
		if (d_rhs[i] != 0.f) {
			if (d_Summ[i] == 0.f)
				d_mlem[i] = d_mlem[i] / d_epps * d_rhs[i];
			else
				d_mlem[i] = d_mlem[i] / d_Summ[i] * d_rhs[i];
		}
		else {
			if (d_Summ[i] != 0.f)
				d_mlem[i] = d_mlem[i] / d_Summ[i] * d_epps;
		}
#endif
		if (d_mlem[i] < d_epps)
			d_mlem[i] = d_epps;
#endif
	}
}

#ifdef PSF
__kernel void Convolution3D(const __global CAST* input, __global CAST* output,
	__constant float* convolution_window, int window_size_x, int window_size_y, int window_size_z) {
	int4 ind = (int4)(get_global_id(0), get_global_id(1), get_global_id(2), 0);
	int4 ind_uus = (int4)(0, 0, 0, 0);
	float result = 0.f;
	const uint Nyx = get_global_size(0) * get_global_size(1);
	//int radius_x = floor((float)window_size_x / 2.0f);
	//int radius_y = floor((float)window_size_y / 2.0f);
	//int radius_z = floor((float)window_size_z / 2.0f);
	int c = 0;

	for (int k = -window_size_z; k <= window_size_z; k++) {
		if (ind.z < window_size_z) {
			if (k < -ind.z)
				ind_uus.z = abs(k) - 1 - ind.z;
			else
				ind_uus.z = k + ind.z;
		}
		else {
			ind_uus.z = ind.z + k;
			if (ind_uus.z >= get_global_size(2))
				ind_uus.z = get_global_size(2) - 1 - (ind_uus.z - get_global_size(2));
			//if (ind_uus.z < 0)
			//	ind_uus.z = ind.z - (k + 1);
		}
		ind_uus.z *= Nyx;
		for (int j = -window_size_y; j <= window_size_y; j++) {
			if (ind.y < window_size_y) {
				if (j < -ind.y)
					ind_uus.y = abs(j) - 1 - ind.y;
				else
					ind_uus.y = j + ind.y;
			}
			else {
				ind_uus.y = ind.y + j;
				if (ind_uus.y >= get_global_size(1))
					ind_uus.y = get_global_size(1) - 1 - (ind_uus.y - get_global_size(1));
				//if (ind_uus.y < 0)
				//	ind_uus.y = ind.y - (j + 1);
			}
			ind_uus.y *= get_global_size(0);
			for (int i = (-window_size_x); i <= window_size_x; i++) {
				//int indx = convert_int(ind.x);
				//indx += i;
				if (ind.x < window_size_x) {
					if (i < -ind.x)
						ind_uus.x = abs(i) - 1 - ind.x;
					else
						ind_uus.x = i + ind.x;
				}
				else {
					ind_uus.x = ind.x + i;
					if (ind_uus.x >= get_global_size(0))
						ind_uus.x = get_global_size(0) - 1 - (ind_uus.x - get_global_size(0));
					//if (ind_uus.x < 0)
					//	ind_uus.x = ind.x - (i + 1);
				}
				//if (indx >= get_global_size(0))
				//	indx = ind.x - i + 1;
				//else if (indx < 0)
				//	indx = abs(convert_int(ind.x)) - 1;
				//int indeksi = indx + ind_uus.y + ind_uus.z;
				int indeksi = ind_uus.x + ind_uus.y + ind_uus.z;
#if defined(ATOMIC) || defined(ATOMIC32)
				float p = convert_float(input[indeksi]) / TH;
#else
				float p = input[indeksi];
#endif
				p *= convolution_window[c];
				result += p;
				c++;
			}
		}
	}
#ifdef ATOMIC
	output[ind.x + ind.y * get_global_size(0) + ind.z * Nyx] = convert_long(result * TH);
#elif defined(ATOMIC32)
	output[ind.x + ind.y * get_global_size(0) + ind.z * Nyx] = convert_int(result * TH);
#else
	output[ind.x + ind.y * get_global_size(0) + ind.z * Nyx] = result;
#endif
}


__kernel void Convolution3D_f(const __global float* input, __global float* output,
	__constant float* convolution_window, int window_size_x, int window_size_y, int window_size_z) {
	int4 ind = (int4)(get_global_id(0), get_global_id(1), get_global_id(2), 0);
	int4 ind_uus = (int4)(0, 0, 0, 0);
	float result = 0.f;
	const uint Nyx = get_global_size(0) * get_global_size(1);
	//int radius_x = floor((float)window_size_x / 2.0f);
	//int radius_y = floor((float)window_size_y / 2.0f);
	//int radius_z = floor((float)window_size_z / 2.0f);
	int c = 0;
	for (int k = -window_size_z; k <= window_size_z; k++) {
		if (ind.z < window_size_z) {
			if (k < -ind.z)
				ind_uus.z = abs(k) - 1 - ind.z;
			else
				ind_uus.z = k + ind.z;
		}
		else {
			ind_uus.z = ind.z + k;
			if (ind_uus.z >= get_global_size(2))
				ind_uus.z = get_global_size(2) - 1 - (ind_uus.z - get_global_size(2));
			//if (ind_uus.z < 0)
			//	ind_uus.z = ind.z - (k + 1);
		}
		ind_uus.z *= Nyx;
		for (int j = -window_size_y; j <= window_size_y; j++) {
			if (ind.y < window_size_y) {
				if (j < -ind.y)
					ind_uus.y = abs(j) - 1 - ind.y;
				else
					ind_uus.y = j + ind.y;
			}
			else {
				ind_uus.y = ind.y + j;
				if (ind_uus.y >= get_global_size(1))
					ind_uus.y = get_global_size(1) - 1 - (ind_uus.y - get_global_size(1));
				//if (ind_uus.y < 0)
				//	ind_uus.y = ind.y - (j + 1);
			}
			ind_uus.y *= get_global_size(0);
			for (int i = (-window_size_x); i <= window_size_x; i++) {
				//int indx = convert_int(ind.x);
				//indx += i;
				if (ind.x < window_size_x) {
					if (i < -ind.x)
						ind_uus.x = abs(i) - 1 - ind.x;
					else
						ind_uus.x = i + ind.x;
				}
				else {
					ind_uus.x = ind.x + i;
					if (ind_uus.x >= get_global_size(0))
						ind_uus.x = get_global_size(0) - 1 - (ind_uus.x - get_global_size(0));
					//if (ind_uus.x < 0)
					//	ind_uus.x = ind.x - (i + 1);
				}
				//if (indx >= get_global_size(0))
				//	indx = ind.x - i + 1;
				//else if (indx < 0)
				//	indx = abs(convert_int(ind.x)) - 1;
				//int indeksi = indx + ind_uus.y + ind_uus.z;
				int indeksi = ind_uus.x + ind_uus.y + ind_uus.z;
				float p = input[indeksi];
				p *= convolution_window[c];
				result += p;
				c++;
			}
		}
	}
	output[ind.x + ind.y * get_global_size(0) + ind.z * Nyx] = result;
}

__kernel void vectorDiv(const __global float* input, __global float* output, const float epps) {
	uint id = get_global_id(0);
	output[id] = output[id] / (input[id] + epps);
}

__kernel void vectorMult(const __global float* input, __global float* output) {
	uint id = get_global_id(0);
	output[id] *= input[id];
}
#endif
#endif

#ifdef NLM_
__kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, 1)))
__kernel void NLM(__global float* grad, const __global float* u, const __global float* u_ref, __constant float* gaussian, const int search_window_x, const int search_window_y,
	const int search_window_z, const int patch_window_x, const int patch_window_y, const int patch_window_z, const uint Nx, const uint Ny, const uint Nz,
	const float h, const float epps, const int Nxy, const int min_x, const int max_x, const int min_y, const int max_y, const int min_z,
	const int max_z, const int type) {

	int x = get_global_id(0);
	int y = get_global_id(1);
	int z = get_global_id(2);
	if (z < min_z || z >= max_z || x < min_x || x >= max_x || y < min_y || y >= max_y)
		return;
	uint n = get_global_id(0) + get_global_id(1) * Nx + get_global_id(2) * Nx * Ny;
	float weight_sum = 0.f;
	float output = 0.f;
	const float uj = u[n];
	for (int k = -search_window_z; k <= search_window_z; k++) {
		const int z_n = z + k;
		for (int j = -search_window_y; j <= search_window_y; j++) {
			const int y_n = y + j;
			for (int i = -search_window_x; i <= search_window_x; i++) {
				const int x_n = x + i;
				const int dim_n = z_n * Nxy + y_n * convert_int(Nx) + x_n;
				const float uk = u[dim_n];
				float distance = 0.f;
				float weight = 0.f;

				for (int pz = -patch_window_z; pz <= patch_window_z; pz++) {
					const int z_k = (z_n + pz) * Nxy;
					const int z_j = (z + pz) * Nxy;
					for (int py = -patch_window_y; py <= patch_window_y; py++) {
						const int y_k = (y_n + py) * convert_int(Nx);
						const int y_j = (y + py) * convert_int(Nx);
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
				weight = exp(-distance / h);
				weight_sum += weight;
				if (type == 2)
					output += weight * uk;
				else if (type == 0) {
					output += (weight * (uj - uk));
				}
				else {
					output += ((weight * (uj - uk)) / sqrt(weight * (uj - uk) * (uj - uk) + epps));
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
__kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, 1)))
__kernel void medianFilter3D(const __global float* grad, __global float* output, const uint Nx, const uint Ny, const uint Nz) {
	int xid = get_global_id(0);
	int yid = get_global_id(1);
	int zid = get_global_id(2);
	if (xid < SEARCH_WINDOW_X || xid >= Nx + SEARCH_WINDOW_X || yid < SEARCH_WINDOW_Y || yid >= Ny + SEARCH_WINDOW_Y || zid < SEARCH_WINDOW_Z || zid >= Nz + SEARCH_WINDOW_Z)
		return;
	int koko = (SEARCH_WINDOW_X * 2 + 1) * (SEARCH_WINDOW_Y * 2 + 1) * (SEARCH_WINDOW_Z * 2 + 1);
	float median[(SEARCH_WINDOW_X * 2 + 1) * (SEARCH_WINDOW_Y * 2 + 1) * (SEARCH_WINDOW_Z * 2 + 1)];
	float medianF[(SEARCH_WINDOW_X * 2 + 1) * (SEARCH_WINDOW_Y * 2 + 1) * (SEARCH_WINDOW_Z * 2 + 1)];
	for (int ll = 0; ll < koko; ll++) {
		medianF[ll] = 0.f;
		median[ll] = 0.f;
	}
	int uu = 0;
	for (int x = -SEARCH_WINDOW_X; x <= SEARCH_WINDOW_X; x++) {
		for (int y = -SEARCH_WINDOW_Y; y <= SEARCH_WINDOW_Y; y++) {
			for (int z = -SEARCH_WINDOW_Z; z <= SEARCH_WINDOW_Z; z++) {
				int pikseli = (xid + x) + (yid + y) * get_global_size(0) + (zid + z) * get_global_size(0) * get_global_size(1);
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
	output[xid + yid * get_global_size(0) + zid * get_global_size(0) * get_global_size(1)] = medianF[koko / 2];
}
#endif