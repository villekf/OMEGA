
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
__constant sampler_t samplerNLM = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_CLAMP;

__kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, 1)))
__kernel void NLM(__global float* restrict grad, __read_only image3d_t restrict u, __constant float* gaussian, const int3 search_window,
	const int3 patch_window, const uint3 N, const float h, const float epps, const int type
#ifdef NLMREF
	, __read_only image3d_t u_ref
#endif
) {

	int3 xyz = { get_global_id(0) , get_global_id(1), get_global_id(2) };
	if (any(xyz < N) || any(xyz >= N))
		return;
	const uint n = xyz.x + xyz.y * N.x + xyz.z * N.x * N.y;
	int4 xyzN = { 0, 0, 0, 0 };
	int4 xyzK = { 0, 0, 0, 0 };
	int4 xyzJ = { 0, 0, 0, 0 };
	float weight_sum = 0.f;
	float output = 0.f;
	const float uj = read_imagef(u, samplerNLM, xyz).x;
	for (int k = -search_window.z; k <= search_window.z; k++) {
		xyzN.z = xyz.z + k;
		for (int j = -search_window.y; j <= search_window.y; j++) {
			xyzN.y = xyz.y + j;
			for (int i = -search_window.x; i <= search_window.x; i++) {
				xyzN.x = xyz.x + i;
				//const int dim_n = z_n * Nxy + y_n * convert_int(N.x) + x_n;
				const float uk = read_imagef(u, samplerNLM, xyzN).x;
				float distance = 0.f;
				float weight = 0.f;

				for (int pz = -patch_window.z; pz <= patch_window.z; pz++) {
					xyzK.z = (xyzN.z + pz);
					xyzJ.z = (xyz.z + pz);
					for (int py = -patch_window.y; py <= patch_window.y; py++) {
						xyzK.y = (xyzN.y + py);
						xyzJ.y = (xyz.y + py);
						int dim_g = (pz + patch_window.z) * (patch_window.x * 2 + 1) * (patch_window.y * 2 + 1) + (py + patch_window.y) * (patch_window.x * 2 + 1);
						for (int px = -patch_window.x; px <= patch_window.x; px++) {
							const float gg = gaussian[dim_g++];
							//const float gg = 1.;
							xyzK.x = (xyzN.x + px);
							xyzJ.x = (xyz.x + px);
#ifdef NLMREF
							const float Pj = read_imagef(u_ref, samplerNLM, xyzK).x;
							const float Pk = read_imagef(u_ref, samplerNLM, xyzJ).x;
#else
							const float Pj = read_imagef(u, samplerNLM, xyzK).x;
							const float Pk = read_imagef(u, samplerNLM, xyzJ).x;
#endif
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

#ifdef RDP
__constant sampler_t samplerRDP = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_CLAMP;

__kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, 1)))
__kernel void RDPKernel(__global float* restrict grad, __read_only image3d_t u, const uint3 N, const float gamma, const float epps) {

	int3 xyz = { get_global_id(0) , get_global_id(1), get_global_id(2) };
	if (any(xyz < N) || any(xyz >= N))
		return;
	const uint n = xyz.x + xyz.y * N.x + xyz.z * N.x * N.y;
	const float uj = read_imagef(u, samplerRDP, xyz).x;
	const float2 ux = { read_imagef(u, samplerRDP, (int4)(xyz.x + 1, xyz.y, xyz.z, 0)).x, read_imagef(u, samplerRDP, (int4)(xyz.x - 1, xyz.y, xyz.z, 0)).x };
	const float2 uy = { read_imagef(u, samplerRDP, (int4)(xyz.x, xyz.y + 1, xyz.z, 0)).x, read_imagef(u, samplerRDP, (int4)(xyz.x, xyz.y - 1, xyz.z, 0)).x };
	const float2 uz = { read_imagef(u, samplerRDP, (int4)(xyz.x, xyz.y, xyz.z + 1, 0)).x, read_imagef(u, samplerRDP, (int4)(xyz.x, xyz.y, xyz.z - 1, 0)).x };
	const float2 uj_ux = uj - ux;
	const float2 uj_uy = uj - uy;
	const float2 uj_uz = uj - uz;

	const float2 divPow2X = (uj + ux + gamma * fabs(uj_ux) + epps);
	const float2 divPow2Y = (uj + uy + gamma * fabs(uj_uy) + epps);
	const float2 divPow2Z = (uj + uz + gamma * fabs(uj_uz) + epps);
	const float2 output = uj_ux * (gamma * fabs(uj_ux) + uj + 3.f * ux + 2.f * epps) / (divPow2X * divPow2X) 
		+ uj_uy * (gamma * fabs(uj_uy) + uj + 3.f * uy + 2.f * epps) / (divPow2Y * divPow2Y)
		+ uj_uz * (gamma * fabs(uj_uz) + uj + 3.f * uz + 2.f * epps) / (divPow2Z * divPow2Z);
	grad[n] = output.x + output.y;
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

#ifdef CPTV
#ifndef DIFFTYPE 
#define DIFFTYPE 0
#endif
//__kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, 1)))
__kernel void CPTVq(__global float* input, const float alpha) {
	ulong idx = get_global_id(0);
	input[idx] = input[idx] * alpha / fmax(fabs(input[idx]), alpha);
}

__kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, 1)))
__kernel void CPTVDivergence(const ulong3 N, const __global float* restrict im, __global float* input) {
	ulong3 xyz = { get_global_id(0) , get_global_id(1), get_global_id(2) };
	if (any(xyz < N) || any(xyz >= N))
		return;
	const ulong x = xyz.x + xyz.y * N.x + xyz.z * N.x * N.y;
	const ulong imDim = (N.x * N.y * N.z);
	ulong xx = x;
	float apuVal = 0.f;
#if DIFFTYPE == 0
	//if (any(xyz == 0) || any(xyz == N - 1)) {
		ulong xh = ((xyz.x - 1) + xyz.y * N.x + xyz.z * N.x * N.y);
		if (xyz.x == 0)
			apuVal += im[xx];
		else
			apuVal += (im[xx] - im[xh]);
		xx *= imDim;
		xh = ((xyz.x) + (xyz.y - 1) * N.x + xyz.z * N.x * N.y) * imDim;
		if (xyz.y == 0)
			apuVal += im[xx];
		else
			apuVal += (im[xx] - im[xh]);
		xx *= imDim;
		xh = ((xyz.x) + (xyz.y) * N.x + (xyz.z - 1) * N.x * N.y) * imDim * 2;
		if (xyz.z == 0)
			apuVal += im[xx];
		else
			apuVal += (im[xx] - im[xh]);
	//}
	//else {
	//	ulong xh = ((xyz.x - 1) + xyz.y * N.x + xyz.z * N.x * N.y);
	//	apuVal -= (im[xx] - im[xh]);
	//	xx *= imDim;
	//	xh *= imDim;
	//	apuVal -= (im[xx] - im[xh]);
	//	xx *= imDim;
	//	xh *= imDim;
	//	apuVal -= (im[xx] - im[xh]);
	//}
#elif DIFFTYPE == 1
#else
#endif
	input[x] -= apuVal;
}

__kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, 1)))
__kernel void CPTVGradient(const ulong3 N, const __global float* restrict im, __global float* input, const float sigma2) {
	ulong3 xyz = { get_global_id(0) , get_global_id(1), get_global_id(2) };
	if (any(xyz < N) || any(xyz >= N))
		return;
	const ulong x = xyz.x + xyz.y * N.x + xyz.z * N.x * N.y;
	const ulong imDim = (N.x * N.y * N.z);
	ulong xx = x;
	float apuVal = 0.f;
#if DIFFTYPE == 0
	float imApu = im[x];
	ulong xh = ((xyz.x + 1) + xyz.y * N.x + xyz.z * N.x * N.y);
	if (xyz.x == N.x - 1)
		apuVal -= imApu;
	else
		apuVal += (im[xh] - imApu);
	input[xx] += apuVal * sigma2;
	apuVal = 0.f;
	xx *= imDim;
	xh = ((xyz.x) + (xyz.y + 1) * N.x + xyz.z * N.x * N.y);
	if (xyz.y == 0)
		apuVal -= imApu;
	else
		apuVal += (im[xh] - imApu);
	input[xx] += apuVal * sigma2;
	apuVal = 0.f;
	xx *= imDim;
	xh = ((xyz.x) + (xyz.y) * N.x + (xyz.z + 1) * N.x * N.y);
	if (xyz.z == 0)
		apuVal -= imApu;
	else
		apuVal += (im[xh] - imApu);
	input[xx] += apuVal * sigma2;
#elif DIFFTYPE == 1
#else
#endif
}
#endif