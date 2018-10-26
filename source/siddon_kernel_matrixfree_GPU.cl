inline void atomicAdd_g_f(volatile __global float *addr, float val) {
	union {
		unsigned int u32;
		float        f32;
	} next, expected, current;
	current.f32    = *addr;
	do {
		expected.f32 = current.f32;
		next.f32     = expected.f32 + val;
		current.u32  = atomic_cmpxchg( (volatile __global unsigned int *)addr, expected.u32, next.u32);
	} while( current.u32 != expected.u32 );
}


__kernel void siddon(__global float* d_rhs, __global float* d_Summ, __global ushort* d_lor, __constant int* d_Nx, __constant int* d_Ny, __constant int* d_Nz,
	__constant float* d_dz, __constant float* d_d, __constant float* d_bz, __constant float* d_bx, __constant float* d_by, __constant float* d_bxf,
	__constant float* d_byf, __constant float* d_bzf, __constant float* d_bxb, __constant float* d_byb, __constant float* d_bzb, __constant float* d_maxxx,
	__constant float* d_maxyy, __constant float* d_minxx, __constant float* d_minyy, __constant float* d_zmax, __constant float* d_NSlices, __global float* d_x,
	__global float* d_y, __global float* d_zdet, __constant int* d_size_x, __constant float* d_xxv, __constant float* d_yyv,
	__constant float* d_xxa, __constant float* d_yya, __global int* d_xyindex, __global int* d_zindex, __constant int* d_TotSinos, __constant int* d_attenuation_correction, 
	__global float* d_atten,  __global float* d_osem, __global float* d_Sino, __constant float* d_epps, __constant int* d_N) {
	int idx = get_global_id(0);
	float xs, xd, ys, yd, zs, zd;
	if (d_xyindex[idx] >= *d_size_x) {
		xs = d_x[d_xyindex[idx]];
		xd = d_x[d_xyindex[idx] - *d_size_x];
		ys = d_y[d_xyindex[idx]];
		yd = d_y[d_xyindex[idx] - *d_size_x];
	}
	else {
		xs = d_x[d_xyindex[idx]];
		xd = d_x[d_xyindex[idx] + *d_size_x];
		ys = d_y[d_xyindex[idx]];
		yd = d_y[d_xyindex[idx] + *d_size_x];
	}
	if (d_zindex[idx] >= *d_TotSinos) {
		zs = d_zdet[d_zindex[idx]];
		zd = d_zdet[d_zindex[idx] - *d_TotSinos];
	}
	else {
		zs = d_zdet[d_zindex[idx]];
		zd = d_zdet[d_zindex[idx] + *d_TotSinos];
	}
	float y_diff = (yd - ys);
	float x_diff = (xd - xs);
	float z_diff = (zd - zs);
	int Np = convert_int(d_lor[idx]);
	int local_ind;
	float local_ele;
	float jelppi = 0.f;
	if (fabs(z_diff) < 1e-8f) {
		int z_loop = convert_int((zs / *d_zmax)*(*d_NSlices - 1.f));
		if (fabs(y_diff) < 1e-8f) {
			if (yd <= *d_maxyy && yd >= *d_minyy) {
				float minvalue = *d_maxyy * 1e5f;
				int apu;
				float start = *d_yya;
				for (int ii = 0; ii < *d_Ny; ii++) {
					float temp = fabs(start - yd);
					start = start + *d_yyv;
					if (temp < minvalue) {
						minvalue = temp;
						apu = ii;
					}
				}
				float templ_ijk = *d_d;
				int tempk = apu * *d_Ny + z_loop * *d_Nx * *d_Ny;
				float temp = *d_d * convert_float(*d_Nx);
				temp = 1.f / temp;
				if (*d_attenuation_correction == 1) {
					for (int iii = 0; iii < Np; iii++) {
						jelppi += (templ_ijk * (-d_atten[tempk + iii]));
					}
					temp = exp(jelppi) * temp;
				}
				templ_ijk = templ_ijk * temp;
				if (d_Sino[idx] != 0.f) {
					float ax = 0.f;
					for (int ii = 0; ii < Np; ii++) {
						local_ind = tempk + ii;
						local_ele = templ_ijk;
						ax += (local_ele * d_osem[local_ind]);
					}
					ax += *d_epps;
					float yax = d_Sino[idx] / ax;
					for (int ii = 0; ii < Np; ii++) {
						local_ind = tempk + ii;
						local_ele = templ_ijk;
						atomicAdd_g_f(&d_rhs[local_ind], (local_ele * yax));
						atomicAdd_g_f(&d_Summ[local_ind], local_ele);
					}
				}
				else {
					for (int ii = 0; ii < Np; ii++) {
						local_ind = tempk + ii;
						local_ele = templ_ijk;
						atomicAdd_g_f(&d_Summ[local_ind], local_ele);
					}
				}
			}
		}
		else if (fabs(x_diff) < 1e-8f) {
			if (xd <= *d_maxxx && xd >= *d_minxx) {
				float minvalue = *d_maxxx * 1e5f;
				int apu;
				float start = *d_xxa;
				for (int ii = 0; ii < *d_Nx; ii++) {
					float temp = fabs(start - xd);
					start = start + *d_xxv;
					if (temp < minvalue) {
						minvalue = temp;
						apu = ii;
					}
				}
				float templ_ijk = *d_d;
				int tempk = apu + z_loop * *d_Nx * *d_Ny;
				float temp = *d_d * convert_float(*d_Ny);
				temp = 1.f / temp;
				if (*d_attenuation_correction == 1) {
					for (int iii = 0; iii < Np; iii++) {
						jelppi += (templ_ijk * (-d_atten[tempk + iii]));
					}
					temp = exp(jelppi) * temp;
				}
				templ_ijk = templ_ijk * temp;
				if (d_Sino[idx] != 0.f) {
					float ax = 0.f;
					for (int ii = 0; ii < Np; ii++) {
						local_ind = tempk + ii;
						local_ele = templ_ijk;
						ax += (local_ele * d_osem[local_ind]);
					}
					ax += *d_epps;
					float yax = d_Sino[idx] / ax;
					for (int ii = 0; ii < Np; ii++) {
						local_ind = tempk + ii;
						local_ele = templ_ijk;
						atomicAdd_g_f(&d_Summ[local_ind], local_ele);
						atomicAdd_g_f(&d_rhs[local_ind], (local_ele * yax));
					}
				}
				else {
					for (int ii = 0; ii < Np; ii++) {
						local_ind = tempk + ii;
						local_ele = templ_ijk;
						atomicAdd_g_f(&d_Summ[local_ind], local_ele);
					}
				}
			}
		}
		else {
			int imin, imax, jmin, jmax, iu, ju, tempi, tempj, tempk;
			float pxt, pyt, tx0, ty0;
			tx0 = (*d_bxf - xs) / (x_diff);
			ty0 = (*d_byf - ys) / (y_diff);
			float txback = (*d_bxb - xs) / (x_diff);
			float tyback = (*d_byb - ys) / (y_diff);
			float txmin = fmin(tx0, txback);
			float txmax = fmax(tx0, txback);
			float tymin = fmin(ty0, tyback);
			float tymax = fmax(ty0, tyback);
			float tmin = fmax(txmin, tymin);
			float tmax = fmin(txmax, tymax);
			if (xs < xd) {
				if (tmin == txmin)
					imin = 1;
				else {
					pxt = xs + tmin*(x_diff);
					imin = convert_int_rtp((pxt - *d_bx) / *d_d);
				}
				if (tmax == txmax)
					imax = *d_Nx;
				else {
					pxt = xs + tmax*(x_diff);
					imax = convert_int_rtn((pxt - *d_bx) / *d_d);
				}
				tx0 = (*d_bx + convert_float(imin) * *d_d - xs) / (x_diff);
				iu = 1;
			}
			else {
				if (tmin == txmin)
					imax = *d_Nx - 1;
				else {
					pxt = xs + tmin*(x_diff);
					imax = convert_int_rtn((pxt - *d_bx) / *d_d);
				}
				if (tmax == txmax)
					imin = 0;
				else {
					pxt = xs + tmax*(x_diff);
					imin = convert_int_rtp((pxt - *d_bx) / *d_d);
				}
				tx0 = (*d_bx + convert_float(imax) * *d_d - xs) / (x_diff);
				iu = -1;
			}
			if (ys < yd) {
				if (tmin == tymin)
					jmin = 1;
				else {
					pyt = ys + tmin*(y_diff);
					jmin = convert_int_rtp((pyt - *d_by) / *d_d);
				}
				if (tmax == tymax)
					jmax = *d_Ny;
				else {
					pyt = ys + tmax*(y_diff);
					jmax = convert_int_rtn((pyt - *d_by) / *d_d);
				}
				ty0 = (*d_by + convert_float(jmin) * *d_d - ys) / (y_diff);
				ju = 1;
			}
			else {
				if (tmin == tymin)
					jmax = *d_Ny - 1;
				else {
					pyt = ys + tmin*(y_diff);
					jmax = convert_int_rtn((pyt - *d_by) / *d_d);
				}
				if (tmax == tymax)
					jmin = 0;
				else {
					pyt = ys + tmax*(y_diff);
					jmin = convert_int_rtp((pyt - *d_by) / *d_d);
				}
				ty0 = (*d_by + convert_float(jmax) * *d_d - ys) / (y_diff);
				ju = -1;
			}
			float pt = ((fmin(tx0, ty0) + tmin) / 2.f);
			tempi = convert_int_rtn(((xs + pt * x_diff) - *d_bx) / *d_d);
			tempj = convert_int_rtn(((ys + pt * y_diff) - *d_by) / *d_d);
			float txu = *d_d / fabs(x_diff);
			float tyu = *d_d / fabs(y_diff);
			float L = sqrt(x_diff*x_diff + y_diff*y_diff);
			tempk = z_loop * *d_Ny * *d_Nx;
			if (ty0 < tx0) {
				if (tmin == ty0)
					tmin -= 1e-7f;
			}
			else {
				if (tmin == tx0)
					tmin -= 1e-7f;
			}
			float tc = tmin;
			float temp = 0.f;
			float tx0_a = tx0, ty0_a = ty0;
			int tempi_a = tempi, tempj_a = tempj, tempk_a = tempk;
			for (int ii = 0; ii < Np; ii++) {
				local_ind = tempj* *d_Nx + tempi + tempk;
				if (local_ind >= *d_N) {
					if (local_ind - *d_Nx >= *d_N)
						local_ind -= (*d_Nx * *d_Ny);
					else if (local_ind - 1 >= *d_N)
						local_ind -= *d_Nx;
					else
						local_ind--;
				}
				else if (local_ind < 0) {
					if (local_ind + *d_Nx < 0)
						local_ind += (*d_Nx * *d_Ny);
					else if (local_ind + 1 < 0)
						local_ind += *d_Nx;
					else
						local_ind++;
				}
				if (tx0 < ty0) {
					local_ele = (tx0 - tc) * L;
					tempi += iu;
					tc = tx0;
					tx0 += txu;
					temp += local_ele;
				}
				else {
					local_ele = (ty0 - tc) * L;
					tempj += ju;
					tc = ty0;
					ty0 += tyu;
					temp += local_ele;
				}
				if (*d_attenuation_correction) {
					jelppi += (local_ele * -d_atten[local_ind]);
				}
			}
			temp = 1.f / temp;
			if (*d_attenuation_correction)
				temp = exp(jelppi) * temp;
			tc = tmin;
			tx0 = tx0_a, ty0 = ty0_a;
			tempi = tempi_a, tempj = tempj_a, tempk = tempk_a;
			if (d_Sino[idx] != 0.f) {
				float ax = 0.f;
				for (int ii = 0; ii < Np; ii++) {
					local_ind = tempj * *d_Nx + tempi + tempk;
					if (local_ind >= *d_N) {
						if (local_ind - *d_Nx >= *d_N)
							local_ind -= (*d_Nx * *d_Ny);
						else if (local_ind - 1 >= *d_N)
							local_ind -= *d_Nx;
						else
							local_ind--;
					}
					else if (local_ind < 0) {
						if (local_ind + *d_Nx < 0)
							local_ind += (*d_Nx * *d_Ny);
						else if (local_ind + 1 < 0)
							local_ind += *d_Nx;
						else
							local_ind++;
					}
					if (tx0 < ty0) {
						local_ele = (tx0 - tc) * L * temp;
						tempi += iu;
						tc = tx0;
						tx0 += txu;
					}
					else {
						local_ele = (ty0 - tc) * L * temp;
						tempj += ju;
						tc = ty0;
						ty0 += tyu;
					}
					ax += (local_ele * d_osem[local_ind]);
				}
				tc = tmin;
				tx0 = tx0_a, ty0 = ty0_a;
				tempi = tempi_a, tempj = tempj_a, tempk = tempk_a;
				ax += *d_epps;
				float yax = d_Sino[idx] / ax;
				for (int ii = 0; ii < Np; ii++) {
					local_ind = tempj * *d_Nx + tempi + tempk;
					if (local_ind >= *d_N) {
						if (local_ind - *d_Nx >= *d_N)
							local_ind -= (*d_Nx * *d_Ny);
						else if (local_ind - 1 >= *d_N)
							local_ind -= *d_Nx;
						else
							local_ind--;
					}
					else if (local_ind < 0) {
						if (local_ind + *d_Nx < 0)
							local_ind += (*d_Nx * *d_Ny);
						else if (local_ind + 1 < 0)
							local_ind += *d_Nx;
						else
							local_ind++;
					}
					if (tx0 < ty0) {
						local_ele = (tx0 - tc) * L * temp;
						tempi += iu;
						tc = tx0;
						tx0 += txu;
					}
					else {
						local_ele = (ty0 - tc) * L * temp;
						tempj += ju;
						tc = ty0;
						ty0 += tyu;
					}
					atomicAdd_g_f(&d_Summ[local_ind], local_ele);
					atomicAdd_g_f(&d_rhs[local_ind], (local_ele * yax));
				}
			}
			else {
				for (int ii = 0; ii < Np; ii++) {
					local_ind = tempj * *d_Nx + tempi + tempk;
					if (local_ind >= *d_N) {
						if (local_ind - *d_Nx >= *d_N)
							local_ind -= (*d_Nx * *d_Ny);
						else if (local_ind - 1 >= *d_N)
							local_ind -= *d_Nx;
						else
							local_ind--;
					}
					else if (local_ind < 0) {
						if (local_ind + *d_Nx < 0)
							local_ind += (*d_Nx * *d_Ny);
						else if (local_ind + 1 < 0)
							local_ind += *d_Nx;
						else
							local_ind++;
					}
					if (tx0 < ty0) {
						local_ele = (tx0 - tc) * L * temp;
						tempi += iu;
						tc = tx0;
						tx0 += txu;
					}
					else {
						local_ele = (ty0 - tc) * L * temp;
						tempj += ju;
						tc = ty0;
						ty0 += tyu;
					}
					atomicAdd_g_f(&d_Summ[local_ind], local_ele);
				}
			}
		}
	}
	else {
		if (fabs(y_diff) < 1e-8f) {
			if (yd <= *d_maxyy && yd >= *d_minyy) {
				float tx0 = (*d_bxf - xs) / (x_diff);
				float tz0 = (*d_bzf - zs) / (z_diff);
				float txback = (*d_bxb - xs) / (x_diff);
				float tzback = (*d_bzb - zs) / (z_diff);
				float txmin = fmin(tx0, txback);
				float txmax = fmax(tx0, txback);
				float tzmin = fmin(tz0, tzback);
				float tzmax = fmax(tz0, tzback);
				float tmin = fmax(txmin, tzmin);
				float tmax = fmin(txmax, tzmax);
				int imin, imax, kmin, kmax;
				float pxt, pzt;
				int iu, ku;
				if (xs < xd) {
					if (tmin == txmin)
						imin = 1;
					else {
						pxt = xs + tmin*(x_diff);
						imin = convert_int_rtp((pxt - *d_bx) / *d_d);
					}
					if (tmax == txmax)
						imax = *d_Nx;
					else {
						pxt = xs + tmax*(x_diff);
						imax = convert_int_rtn((pxt - *d_bx) / *d_d);
					}
					tx0 = (*d_bx + convert_float(imin) * *d_d - xs) / (x_diff);
					iu = 1;
				}
				else if (xs > xd) {
					if (tmin == txmin)
						imax = *d_Nx - 1;
					else {
						pxt = xs + tmin*(x_diff);
						imax = convert_int_rtn((pxt - *d_bx) / *d_d);
					}
					if (tmax == txmax)
						imin = 0;
					else {
						pxt = xs + tmax*(x_diff);
						imin = convert_int_rtp((pxt - *d_bx) / *d_d);
					}
						tx0 = (*d_bx + convert_float(imax) * *d_d - xs) / (x_diff);
					iu = -1;
				}
				if (zs < zd) {
					if (tmin == tzmin)
						kmin = 1;
						else {
						pzt = zs + tmin*(z_diff);
						kmin = convert_int_rtp((pzt - *d_bz) / *d_dz);
					}
					if (tmax == tzmax)
						kmax = *d_Nz;
					else {
						pzt = zs + tmax*(z_diff);
						kmax = convert_int_rtn((pzt - *d_bz) / *d_dz);
					}
					tz0 = (*d_bz + convert_float(kmin) * *d_dz - zs) / (z_diff);
					ku = 1;
				}
				else if (zs > zd) {
					if (tmin == tzmin)
						kmax = *d_Nz - 1;
					else {
						pzt = zs + tmin*(z_diff);
						kmax = convert_int_rtn((pzt - *d_bz) / *d_dz);
					}
					if (tmax == tzmax)
						kmin = 0;
					else {
						pzt = zs + tmax*(z_diff);
						kmin = convert_int_rtp((pzt - *d_bz) / *d_dz);
					}
					tz0 = (*d_bz + convert_float(kmax) * *d_dz - zs) / (z_diff);
					ku = -1;
				}
				float L = sqrt((x_diff*x_diff + z_diff*z_diff));
				int	tempi, tempj, tempk;
				float apu2, apu1;
				float pt = ((fmin(tx0, tz0) + tmin) / 2.f);
				tempi = convert_int_rtn(((xs + pt * x_diff) - *d_bx) / *d_d);
				tempk = convert_int_rtn(((zs + pt * z_diff) - *d_bz) / *d_dz);
				float txu = *d_d / fabs(x_diff);
				float tzu = *d_dz / fabs(z_diff);
				float start = *d_yya;
				for (int ii = 0; ii < Np; ii++) {
					apu1 = fabs(start - yd);
					start = start + *d_yyv;
					if ((ii > 0 && apu1 < apu2) || ii == 0) {
						tempj = ii;
						apu2 = apu1;
					}
				}
				tempj = tempj * *d_Nx;
				if (tz0 < tx0) {
					if (tmin == tz0)
						tmin -= 1e-7f;
				}
				else {
					if (tmin == tx0)
						tmin -= 1e-7f;
				}
				float tc = tmin;
				float temp = 0.f;
				float tx0_a = tx0, tz0_a = tz0;
				int tempi_a = tempi, tempj_a = tempj, tempk_a = tempk;
				for (int ii = 0; ii < Np; ii++) {
					local_ind = tempj + tempi + *d_Nx * *d_Ny * tempk;
					if (local_ind >= *d_N) {
						if (local_ind - *d_Nx >= *d_N)
							local_ind -= (*d_Nx - *d_Ny);
						else if (local_ind - 1 >= *d_N)
							local_ind -= *d_Nx;
						else
							local_ind--;
					}
					else if (local_ind < 0) {
						if (local_ind + *d_Nx < 0)
							local_ind += (*d_Nx * *d_Ny);
						else if (local_ind + 1 < 0)
							local_ind += *d_Nx;
						else
							local_ind++;
					}
					if (tx0 < tz0) {
						local_ele = (tx0 - tc) * L;
						tempi += iu;
						tc = tx0;
						tx0 += txu;
						temp += local_ele;
					}
					else {
						local_ele = (tz0 - tc) * L;
						tempk += ku;
						tc = tz0;
						tz0 += tzu;
						temp += local_ele;
					}
					if (*d_attenuation_correction) {
						jelppi += (local_ele * -d_atten[local_ind]);
					}
				}
				temp = 1.f / temp;
				if (*d_attenuation_correction)
					temp = exp(jelppi) * temp;
				tc = tmin;
				tx0 = tx0_a, tz0 = tz0_a;
				tempi = tempi_a, tempj = tempj_a, tempk = tempk_a;
				if (d_Sino[idx] != 0.f) {
					float ax = 0.f;
					for (int ii = 0; ii < Np; ii++) {
						local_ind = tempj + tempi + *d_Nx * *d_Ny * tempk;
						if (local_ind >= *d_N) {
							if (local_ind - *d_Nx >= *d_N)
								local_ind -= (*d_Nx * *d_Ny);
							else if (local_ind - 1 >= *d_N)
								local_ind -= *d_Nx;
							else
								local_ind--;
						}
						else if (local_ind < 0) {
							if (local_ind + *d_Nx < 0)
								local_ind += (*d_Nx * *d_Ny);
							else if (local_ind + 1 < 0)
								local_ind += *d_Nx;
							else
								local_ind++;
						}
						if (tx0 < tz0) {
							local_ele = (tx0 - tc) * L * temp;
							tempi += iu;
							tc = tx0;
							tx0 += txu;
						}
						else {
							local_ele = (tz0 - tc) * L * temp;
							tempk += ku;
							tc = tz0;
							tz0 += tzu;
						}
						ax += (local_ele * d_osem[local_ind]);
					}
					tc = tmin;
					tx0 = tx0_a, tz0 = tz0_a;
					tempi = tempi_a, tempj = tempj_a, tempk = tempk_a;
					ax += *d_epps;
					float yax = d_Sino[idx] / ax;
					for (int ii = 0; ii < Np; ii++) {
						local_ind = tempj + tempi + *d_Nx * *d_Ny * tempk;
						if (local_ind >= *d_N) {
							if (local_ind - *d_Nx >= *d_N)
								local_ind -= (*d_Nx * *d_Ny);
							else if (local_ind - 1 >= *d_N)
								local_ind -= *d_Nx;
							else
								local_ind--;
						}
						else if (local_ind < 0) {
							if (local_ind + *d_Nx < 0)
								local_ind += (*d_Nx * *d_Ny);
							else if (local_ind + 1 < 0)
								local_ind += *d_Nx;
							else
								local_ind++;
						}
						if (tx0 < tz0) {
							local_ele = (tx0 - tc) * L * temp;
							tempi += iu;
							tc = tx0;
							tx0 += txu;
						}
						else {
							local_ele = (tz0 - tc) * L * temp;
							tempk += ku;
							tc = tz0;
							tz0 += tzu;
						}
						atomicAdd_g_f(&d_Summ[local_ind], local_ele);
						atomicAdd_g_f(&d_rhs[local_ind], (local_ele * yax));
					}
				}
				else {
					for (int ii = 0; ii < Np; ii++) {
						local_ind = tempj + tempi + *d_Nx * *d_Ny * tempk;
						if (local_ind >= *d_N) {
							if (local_ind - *d_Nx >= *d_N)
								local_ind -= (*d_Nx * *d_Ny);
							else if (local_ind - 1 >= *d_N)
								local_ind -= *d_Nx;
							else
								local_ind--;
						}
						else if (local_ind < 0) {
							if (local_ind + *d_Nx < 0)
								local_ind += (*d_Nx * *d_Ny);
							else if (local_ind + 1 < 0)
								local_ind += *d_Nx;
							else
								local_ind++;
						}
						if (tx0 < tz0) {
							local_ele = (tx0 - tc) * L * temp;
							tempi += iu;
							tc = tx0;
							tx0 += txu;
						}
						else {
							local_ele = (tz0 - tc) * L * temp;
							tempk += ku;
							tc = tz0;
							tz0 += tzu;
						}
						atomicAdd_g_f(&d_Summ[local_ind], local_ele);
					}
				}
			}
		}
		else if (fabs(x_diff) < 1e-8f) {
			if (xd <= *d_maxxx && xd >= *d_minxx) {
				float ty0 = (*d_byf - ys) / (y_diff);
				float tyback = (*d_byb - ys) / (y_diff);
				float tz0 = (*d_bzf - zs) / (z_diff);
				float tzback = (*d_bzb - zs) / (z_diff);
				float tzmin = fmin(tz0, tzback);
				float tzmax = fmax(tz0, tzback);
				float tymin = fmin(ty0, tyback);
				float tymax = fmax(ty0, tyback);
				float tmin = fmax(tymin, tzmin);
				float tmax = fmin(tymax, tzmax);
				int jmin, jmax, kmin, kmax, ku, ju;
				float pyt, pzt;
				if (ys < yd) {
					if (tmin == tymin)
						jmin = 1;
					else {
						pyt = ys + tmin*(y_diff);
						jmin = convert_int_rtp((pyt - *d_by) / *d_d);
					}
					if (tmax == tymax)
						jmax = *d_Ny;
					else {
						pyt = ys + tmax*(y_diff);
						jmax = convert_int_rtn((pyt - *d_by) / *d_d);
					}
					ty0 = (*d_by + convert_float(jmin) * *d_d - ys) / (y_diff);
					ju = 1;
				}
				else if (ys > yd) {
					if (tmin == tymin)
						jmax = *d_Ny - 1;
					else {
						pyt = ys + tmin*(y_diff);
						jmax = convert_int_rtn((pyt - *d_by) / *d_d);
					}
					if (tmax == tymax)
						jmin = 0;
					else {
						pyt = ys + tmax*(y_diff);
						jmin = convert_int_rtp((pyt - *d_by) / *d_d);
					}
					ty0 = (*d_by + convert_float(jmax) * *d_d - ys) / (y_diff);
					ju = -1;
				}
				if (zs < zd) {
					if (tmin == tzmin)
						kmin = 1;
					else {
						pzt = zs + tmin*(z_diff);
						kmin = convert_int_rtp((pzt - *d_bz) / *d_dz);
					}
					if (tmax == tzmax)
						kmax = *d_Nz;
					else {
						pzt = zs + tmax*(z_diff);
						kmax = convert_int_rtn((pzt - *d_bz) / *d_dz);
					}
					tz0 = (*d_bz + convert_float(kmin) * *d_dz - zs) / (z_diff);
					ku = 1;
				}
				else if (zs > zd) {
					if (tmin == tzmin)
						kmax = *d_Nz - 1;
					else {
						pzt = zs + tmin*(z_diff);
						kmax = convert_int_rtn((pzt - *d_bz) / *d_dz);
					}
					if (tmax == tzmax)
						kmin = 0;
					else {
						pzt = zs + tmax*(z_diff);
						kmin = convert_int_rtp((pzt - *d_bz) / *d_dz);
					}
					tz0 = (*d_bz + convert_float(kmax) * *d_dz - zs) / (z_diff);
					ku = -1;
				}
				float L = sqrt((y_diff*y_diff + z_diff*z_diff));
				int tempi, tempj, tempk;
				float apu2, apu1, apu;
				float pt = ((min(tz0, ty0) + tmin) / 2.);
				tempk = convert_int_rtn(((zs + pt * z_diff) - *d_bz) / *d_dz);
				tempj = convert_int_rtn(((ys + pt * y_diff) - *d_by) / *d_d);
				float tzu = *d_dz / fabs(z_diff);
				float tyu = *d_d / fabs(y_diff);
				if (tz0 < ty0) {
					if (tmin == tz0)
						tmin -= 1e-7f;
				}
				else {
					if (tmin == ty0)
						tmin -= 1e-7f;
				}
				float tc = tmin;
				float temp = 0.f;
				float start = *d_xxa;
				for (int ii = 0; ii < *d_Nx; ii++) {
					apu1 = fabs(*d_xxa - xd);
					start = start + *d_xxv;
					if ((ii > 0 && apu1 < apu2) || ii == 0) {
						tempi = ii;
						apu2 = apu1;
					}
				}
				float tz0_a = tz0, ty0_a = ty0;
				int tempi_a = tempi, tempj_a = tempj, tempk_a = tempk;
				for (int ii = 0; ii < Np; ii++) {
					local_ind = tempj * *d_Nx + tempi + *d_Nx * *d_Ny * tempk;
					if (local_ind >= *d_N) {
						if (local_ind - *d_Nx >= *d_N)
							local_ind -= (*d_Nx * *d_Ny);
						else if (local_ind - 1 >= *d_N)
							local_ind -= *d_Nx;
						else
							local_ind--;
					}
					else if (local_ind < 0) {
						if (local_ind + *d_Nx < 0)
							local_ind += (*d_Nx * *d_Ny);
						else if (local_ind + 1 < 0)
							local_ind += *d_Nx;
						else
							local_ind++;
					}
					if (tz0 < ty0) {
						local_ele = (tz0 - tc) * L;
						tempk += ku;
						tc = tz0;
						tz0 += tzu;
						temp += local_ele;
					}
					else {
						local_ele = (ty0 - tc) * L;
						tempj += ju;
						tc = ty0;
						ty0 += tyu;
						temp += local_ele;
					}
					if (*d_attenuation_correction) {
						jelppi += (local_ele * -d_atten[local_ind]);
					}
				}
				temp = 1.f / temp;
				if (*d_attenuation_correction)
					temp = exp(jelppi) * temp;
				tc = tmin;
				tz0 = tz0_a, ty0 = ty0_a;
				tempi = tempi_a, tempj = tempj_a, tempk = tempk_a;
				if (d_Sino[idx] != 0.f) {
					float ax = 0.f;
					for (int ii = 0; ii < Np; ii++) {
						local_ind = tempj * *d_Nx + tempi + *d_Nx * *d_Ny * tempk;
						if (local_ind >= *d_N) {
							if (local_ind - *d_Nx >= *d_N)
								local_ind -= (*d_Nx * *d_Ny);
							else if (local_ind - 1 >= *d_N)
								local_ind -= *d_Nx;
							else
								local_ind--;
						}
						else if (local_ind < 0) {
							if (local_ind + *d_Nx < 0)
								local_ind += (*d_Nx * *d_Ny);
							else if (local_ind + 1 < 0)
								local_ind += *d_Nx;
							else
								local_ind++;
						}
						if (tz0 < ty0) {
							local_ele = (tz0 - tc) * L * temp;
							tempk += ku;
							tc = tz0;
							tz0 += tzu;
						}
						else {
							local_ele = (ty0 - tc) * L * temp;
							tempj += ju;
							tc = ty0;
							ty0 += tyu;
						}
						ax += (local_ele * d_osem[local_ind]);
					}
					tc = tmin;
					tz0 = tz0_a, ty0 = ty0_a;
					tempi = tempi_a, tempj = tempj_a, tempk = tempk_a;
					ax += *d_epps;
					float yax = d_Sino[idx] / ax;
					for (int ii = 0; ii < Np; ii++) {
						local_ind = tempj * *d_Nx + tempi + *d_Nx * *d_Ny * tempk;
						if (local_ind >= *d_N) {
							if (local_ind - *d_Nx >= *d_N)
								local_ind -= (*d_Nx * *d_Ny);
							else if (local_ind - 1 >= *d_N)
								local_ind -= *d_Nx;
							else
								local_ind--;
						}
						else if (local_ind < 0) {
							if (local_ind + *d_Nx < 0)
								local_ind += (*d_Nx * *d_Ny);
							else if (local_ind + 1 < 0)
								local_ind += *d_Nx;
							else
								local_ind++;
						}
						if (tz0 < ty0) {
							local_ele = (tz0 - tc) * L * temp;
							tempk += ku;
							tc = tz0;
							tz0 += tzu;
						}
						else {
							local_ele = (ty0 - tc) * L * temp;
							tempj += ju;
							tc = ty0;
							ty0 += tyu;
						}
						atomicAdd_g_f(&d_Summ[local_ind], local_ele);
						atomicAdd_g_f(&d_rhs[local_ind], (local_ele * yax));
					}
				}
				else {
					for (int ii = 0; ii < Np; ii++) {
						local_ind = tempj * *d_Nx + tempi + *d_Nx * *d_Ny * tempk;
						if (local_ind >= *d_N) {
							if (local_ind - *d_Nx >= *d_N)
								local_ind -= (*d_Nx * *d_Ny);
							else if (local_ind - 1 >= *d_N)
								local_ind -= *d_Nx;
							else
								local_ind--;
						}
						else if (local_ind < 0) {
							if (local_ind + *d_Nx < 0)
								local_ind += (*d_Nx * *d_Ny);
							else if (local_ind + 1 < 0)
								local_ind += *d_Nx;
							else
								local_ind++;
						}
						if (tz0 < ty0) {
							local_ele = (tz0 - tc) * L * temp;
							tempk += ku;
							tc = tz0;
							tz0 += tzu;
						}
						else {
							local_ele = (ty0 - tc) * L * temp;
							tempj += ju;
							tc = ty0;
							ty0 += tyu;
						}
						atomicAdd_g_f(&d_Summ[local_ind], local_ele);
					}
				}
			}
		}
		else {
			float tx0 = (*d_bxf - xs) / (x_diff);
			float tz0 = (*d_bzf - zs) / (z_diff);
			float txback = (*d_bxb - xs) / (x_diff);
			float tzback = (*d_bzb - zs) / (z_diff);
			float ty0 = (*d_byf - ys) / (y_diff);
			float tyback = (*d_byb - ys) / (y_diff);
			float txmin = fmin(tx0, txback);
			float txmax = fmax(tx0, txback);
			float tymin = fmin(ty0, tyback);
			float tymax = fmax(ty0, tyback);
			float tzmin = fmin(tz0, tzback);
			float tzmax = fmax(tz0, tzback);
			float tmin = fmax(max(txmin, tzmin), tymin);
			float tmax = fmin(min(txmax, tzmax), tymax);
			int imin, imax, jmin, jmax, kmin, kmax, iu, ju, ku;
			float pxt, pyt, pzt;
			if (xs < xd) {
				if (tmin == txmin)
					imin = 1;
				else {
					pxt = xs + tmin*(x_diff);
					imin = convert_int_rtp((pxt - *d_bx) / *d_d);
				}
				if (tmax == txmax)
					imax = *d_Nx;
				else {
					pxt = xs + tmax*(x_diff);
					imax = convert_int_rtn((pxt - *d_bx) / *d_d);
				}
				tx0 = (*d_bx + convert_float(imin) * *d_d - xs) / (x_diff);
				iu = 1;
			}
			else if (xs > xd) {
				if (tmin == txmin)
					imax = *d_Nx - 1;
				else {
					pxt = xs + tmin*(x_diff);
					imax = convert_int_rtn((pxt - *d_bx) / *d_d);
				}
				if (tmax == txmax)
					imin = 0;
				else {
					pxt = xs + tmax*(x_diff);
					imin = convert_int_rtp((pxt - *d_bx) / *d_d);
				}
				tx0 = (*d_bx + convert_float(imax) * *d_d - xs) / (x_diff);
				iu = -1;
			}
			if (ys < yd) {
				if (tmin == tymin)
					jmin = 1;
				else {
					pyt = ys + tmin*(y_diff);
					jmin = convert_int_rtp((pyt - *d_by) / *d_d);
				}
				if (tmax == tymax)
					jmax = *d_Ny;
				else {
					pyt = ys + tmax*(y_diff);
					jmax = convert_int_rtn((pyt - *d_by) / *d_d);
				}
				ty0 = (*d_by + convert_float(jmin) * *d_d - ys) / (y_diff);
				ju = 1;
			}
			else if (ys > yd) {
				if (tmin == tymin)
					jmax = *d_Ny - 1;
				else {
					pyt = ys + tmin*(y_diff);
					jmax = convert_int_rtn((pyt - *d_by) / *d_d);
				}
				if (tmax == tymax)
					jmin = 0;
				else {
					pyt = ys + tmax*(y_diff);
					jmin = convert_int_rtp((pyt - *d_by) / *d_d);
				}
				ty0 = (*d_by + convert_float(jmax) * *d_d - ys) / (y_diff);
				ju = -1;
			}
			if (zs < zd) {
				if (tmin == tzmin)
					kmin = 1;
				else {
					pzt = zs + tmin*(z_diff);
					kmin = convert_int_rtp((pzt - *d_bz) / *d_dz);
				}
				if (tmax == tzmax)
					kmax = *d_Nz;
				else {
					pzt = zs + tmax*(z_diff);
					kmax = convert_int_rtn((pzt - *d_bz) / *d_dz);
				}
				tz0 = (*d_bz + convert_float(kmin) * *d_dz - zs) / (z_diff);
				ku = 1;
			}
			else if (zs > zd) {
				if (tmin == tzmin)
					kmax = *d_Nz - 1;
				else {
					pzt = zs + tmin*(z_diff);
					kmax = convert_int_rtn((pzt - *d_bz) / *d_dz);
				}
				if (tmax == tzmax)
					kmin = 0;
				else {
					pzt = zs + tmax*(z_diff);
					kmin = convert_int_rtp((pzt - *d_bz) / *d_dz);
				}
				tz0 = (*d_bz + convert_float(kmax) * *d_dz - zs) / (z_diff);
				ku = -1;
			}
			float L = sqrt(x_diff*x_diff + z_diff*z_diff + y_diff*y_diff);
			int tempi, tempj, tempk;
			float pt = ((fmin(fmin(tz0, ty0), tx0) + tmin) / 2.f);
			tempk = convert_int_rtn(((zs + pt * z_diff) - *d_bz) / *d_dz);
			tempj = convert_int_rtn(((ys + pt * y_diff) - *d_by) / *d_d);
			tempi = convert_int_rtn(((xs + pt * x_diff) - *d_bx) / *d_d);
			float tzu = *d_dz / fabs(z_diff);
			float tyu = *d_d / fabs(y_diff);
			float txu = *d_d / fabs(x_diff);
			if (tz0 < ty0 && tz0 < tx0) {
				if (tmin == tz0)
					tmin -= 1e-7f;
			}
			else if (ty0 < tx0) {
				if (tmin == ty0)
					tmin -= 1e-7f;
			}
			else {
				if (tmin == tx0)
					tmin -= 1e-7f;
			}
			float tc = tmin;
			float temp = 0.f;
			float tx0_a = tx0, ty0_a = ty0, tz0_a = tz0;
			int tempi_a = tempi, tempj_a = tempj, tempk_a = tempk;
			for (int ii = 0; ii < Np; ii++) {
				local_ind = tempj * *d_Nx + tempi + *d_Nx * *d_Ny * tempk;
				if (local_ind >= *d_N) {
					if (local_ind - *d_Nx >= *d_N)
						local_ind -= (*d_Nx * *d_Ny);
					else if (local_ind - 1 >= *d_N)
						local_ind -= *d_Nx;
					else
						local_ind--;
				}
				else if (local_ind < 0) {
					if (local_ind + *d_Nx < 0)
						local_ind += (*d_Nx * *d_Ny);
					else if (local_ind + 1 < 0)
						local_ind += *d_Nx;
					else
						local_ind++;
				}
				if (tz0 < ty0 && tz0 < tx0) {
					local_ele = (tz0 - tc) * L;
					tempk += ku;
					tc = tz0;
					tz0 += tzu;
					temp += local_ele;
				}
				else if (ty0 < tx0) {
					local_ele = (ty0 - tc) * L;
					tempj += ju;
					tc = ty0;
					ty0 += tyu;
					temp += local_ele;
				}
				else {
					local_ele = (tx0 - tc) * L;
					tempi += iu;
					tc = tx0;
					tx0 += txu;
					temp += local_ele;
				}
				if (*d_attenuation_correction) {
					jelppi += (local_ele * -d_atten[local_ind]);
				}
			}
			temp = 1.f / temp;
			if (*d_attenuation_correction)
				temp = exp(jelppi) * temp;
			tc = tmin;
			tx0 = tx0_a, ty0 = ty0_a, tz0 = tz0_a;
			tempi = tempi_a, tempj = tempj_a, tempk = tempk_a;
			if (d_Sino[idx] != 0.f) {
				float ax = 0.f;
				for (int ii = 0; ii < Np; ii++) {
					local_ind = tempj * *d_Nx + tempi + *d_Nx * *d_Ny * tempk;
					if (local_ind >= *d_N) {
						if (local_ind - *d_Nx >= *d_N)
							local_ind -= (*d_Nx * *d_Ny);
						else if (local_ind - 1 >= *d_N)
							local_ind -= *d_Nx;
						else
							local_ind--;
					}
					else if (local_ind < 0) {
						if (local_ind + *d_Nx < 0)
							local_ind += (*d_Nx * *d_Ny);
						else if (local_ind + 1 < 0)
							local_ind += *d_Nx;
						else
							local_ind++;
					}
					if (tz0 < ty0 && tz0 < tx0) {
						local_ele = (tz0 - tc) * L * temp;
						tempk += ku;
						tc = tz0;
						tz0 += tzu;
					}
					else if (ty0 < tx0) {
						local_ele = (ty0 - tc) * L * temp;
						tempj += ju;
						tc = ty0;
						ty0 += tyu;
					}
					else {
						local_ele = (tx0 - tc) * L * temp;
						tempi += iu;
						tc = tx0;
						tx0 += txu;
					}
					ax += (local_ele * d_osem[local_ind]);
				}
				tc = tmin;
				tx0 = tx0_a, ty0 = ty0_a, tz0 = tz0_a;
				tempi = tempi_a, tempj = tempj_a, tempk = tempk_a;
				ax += *d_epps;
				float yax = d_Sino[idx] / ax;
				for (int ii = 0; ii < Np; ii++) {
					local_ind = tempj * *d_Nx + tempi + *d_Nx * *d_Ny * tempk;
					if (local_ind >= *d_N) {
						if (local_ind - *d_Nx >= *d_N)
							local_ind -= (*d_Nx * *d_Ny);
						else if (local_ind - 1 >= *d_N)
							local_ind -= *d_Nx;
						else
							local_ind--;
					}
					else if (local_ind < 0) {
						if (local_ind + *d_Nx < 0)
							local_ind += (*d_Nx * *d_Ny);
						else if (local_ind + 1 < 0)
							local_ind += *d_Nx;
						else
							local_ind++;
					}
					if (tz0 < ty0 && tz0 < tx0) {
						local_ele = (tz0 - tc) * L * temp;
						tempk += ku;
						tc = tz0;
						tz0 += tzu;
					}
					else if (ty0 < tx0) {
						local_ele = (ty0 - tc) * L * temp;
						tempj += ju;
						tc = ty0;
						ty0 += tyu;
					}
					else {
						local_ele = (tx0 - tc) * L * temp;
						tempi += iu;
						tc = tx0;
						tx0 += txu;
					}
					//barrier(CLK_GLOBAL_MEM_FENCE);
					atomicAdd_g_f(&d_Summ[local_ind], local_ele);
					atomicAdd_g_f(&d_rhs[local_ind], (local_ele * yax));
				}
			}
			else {
				for (int ii = 0; ii < Np; ii++) {
					local_ind = tempj * *d_Nx + tempi + *d_Nx * *d_Ny * tempk;
					if (local_ind >= *d_N) {
						if (local_ind - *d_Nx >= *d_N)
							local_ind -= (*d_Nx * *d_Ny);
						else if (local_ind - 1 >= *d_N)
							local_ind -= *d_Nx;
						else
							local_ind--;
					}
					else if (local_ind < 0) {
						if (local_ind + *d_Nx < 0)
							local_ind += (*d_Nx * *d_Ny);
						else if (local_ind + 1 < 0)
							local_ind += *d_Nx;
						else
							local_ind++;
					}
					if (tz0 < ty0 && tz0 < tx0) {
						local_ele = (tz0 - tc) * L * temp;
						tempk += ku;
						tc = tz0;
						tz0 += tzu;
					}
					else if (ty0 < tx0) {
						local_ele = (ty0 - tc) * L * temp;
						tempj += ju;
						tc = ty0;
						ty0 += tyu;
					}
					else {
						local_ele = (tx0 - tc) * L * temp;
						tempi += iu;
						tc = tx0;
						tx0 += txu;
					}
					atomicAdd_g_f(&d_Summ[local_ind], local_ele);
				}
			}
		}
	}
}