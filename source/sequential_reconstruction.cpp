#include "mex.h"
#include <vector>
#include <algorithm>
#include <time.h>


using namespace std;

void improved_siddon(const size_t loop_var_par, const int size_x, const double zmax, double* Summ, double* rhs, const double maxyy, const double minyy, 
	const double maxxx, const double minxx, const vector<double>& xx_vec, const vector<double>& iij_vec, const vector<double>& jjk_vec, 
	const vector<double>& kkj_vec, const vector<double>& yy_vec, const double* atten, const double* x, const double* y, const double* z_det, const int NSlices, 
	const int Nx, const int Ny, const int Nz, const double d, const double dz, const double bx, const double by, const double bz, const bool attenuation_correction, 
	const uint16_t* lor1, const uint32_t* xy_index, const uint32_t* z_index, const int TotSinos, const int maksimi, const double epps, const double* Sino, 
	double* osem_apu, const int N) {

	int oo = -1;

	bool pass = false;

	const double bxf = bx + iij_vec.front() * d;
	const double byf = by + jjk_vec.front() * d;
	const double bxb = bx + iij_vec.back() * d;
	const double byb = by + jjk_vec.back() * d;
	const double bzf = bz + kkj_vec.front() * dz;
	const double bzb = bz + kkj_vec.back() * dz;

	vector<uint32_t> indices(maksimi, 0);
	vector<double> elements(maksimi, 0.);

	for (int lo = 1; lo <= loop_var_par; lo++) {
			
		oo++;

		pass = false;

		double xs, xd, ys, yd, zs, zd;
		if (xy_index[oo] >= size_x) {
			xs = x[xy_index[oo]];
			xd = x[xy_index[oo] - size_x];
			ys = y[xy_index[oo]];
			yd = y[xy_index[oo] - size_x];
		}
		else {
			xs = x[xy_index[oo]];
			xd = x[xy_index[oo] + size_x];
			ys = y[xy_index[oo]];
			yd = y[xy_index[oo] + size_x];
		}
		if (z_index[oo] >= TotSinos) {
			zs = z_det[z_index[oo]];
			zd = z_det[z_index[oo] - TotSinos];
		}
		else {
			zs = z_det[z_index[oo]];
			zd = z_det[z_index[oo] + TotSinos];
		}
		

		double y_diff = (yd - ys);
		double x_diff = (xd - xs);
		double z_diff = (zd - zs);

		unsigned int Np = static_cast<int>(lor1[oo]);
		

		if (fabs(z_diff) < 1e-8) {

			int z_loop = static_cast<int>((zs / zmax)*(static_cast<double>(NSlices) - 1.));

			if (fabs(y_diff) < 1e-8) {

				if (yd <= maxyy && yd >= minyy) {
					double minvalue = maxyy * 1e5;
					int apu;
					for (int ii = 0; ii < Ny; ii++) {
						double temp = fabs(yy_vec[ii] - yd);
						if (temp < minvalue) {
							minvalue = temp;
							apu = ii;
						}
					}
					double templ_ijk = d;
					int tempk = apu * Ny + z_loop * Nx * Ny;
					double temp = d * static_cast<double>(Nx);
					temp = 1. / temp;

					if (attenuation_correction) {

						double jelppi = 0.;

						for (unsigned int iii = 0; iii < Np; iii++) {
							jelppi += templ_ijk * -atten[tempk + iii];
						}
						temp = exp(jelppi) * temp;
					}

					double element = (templ_ijk * temp);

					for (unsigned int ii = 0; ii < Np; ii++) {
						indices[ii] = tempk + ii;
						elements[ii] = element;
					}
					pass = true;
				}
			}
			else if (fabs(x_diff) < 1e-8) {

				if (xd <= maxxx && xd >= minxx) {
					double minvalue = maxxx * 1e5;
					int apu;
					for (int ii = 0; ii < Nx; ii++) {
						double temp = fabs(xx_vec[ii] - xd);
						if (temp < minvalue) {
							minvalue = temp;
							apu = ii;
						}
					}
					double templ_ijk = d;
					int tempk = apu + z_loop * Nx * Ny;
					double temp = d * static_cast<double>(Ny);
					temp = 1. / temp;

					if (attenuation_correction) {

						double jelppi = 0.;

						for (unsigned int iii = 0; iii < Np; iii++) {
							jelppi += templ_ijk * -atten[tempk + iii * Ny];
						}
						temp = exp(jelppi) * temp;

					}

					double element = (templ_ijk * temp);

					for (unsigned int ii = 0; ii < Np; ii++) {
						indices[ii] = tempk + ii * Ny;
						elements[ii] = element;
					}
					pass = true;
				}
			}
			else {
				double tx0 = (bxf - xs) / (x_diff);
				double ty0 = (byf - ys) / (y_diff);
				double txback = (bxb - xs) / (x_diff);
				double tyback = (byb - ys) / (y_diff);

				double txmin = min(tx0, txback);
				double txmax = max(tx0, txback);
				double tymin = min(ty0, tyback);
				double tymax = max(ty0, tyback);

				double tmin = max(txmin, tymin);
				double tmax = min(txmax, tymax);

				int imin, imax, jmin, jmax;
				double pxt, pyt;

				int iu, ju;

				if (xs < xd) {
					if (tmin == txmin)
						imin = 1;
					else {
						pxt = xs + tmin*(x_diff);
						imin = static_cast<int>(ceil((pxt - bx) / d));
					}
					if (tmax == txmax)
						imax = Nx;
					else {
						pxt = xs + tmax*(x_diff);
						imax = static_cast<int>(floor((pxt - bx) / d));
					}
					tx0 = (bx + static_cast<double>(imin) * d - xs) / (x_diff);
					iu = 1;
				}
				else if (xs > xd) {
					if (tmin == txmin)
						imax = Nx - 1;
					else {
						pxt = xs + tmin*(x_diff);
						imax = static_cast<int>(floor((pxt - bx) / d));
					}
					if (tmax == txmax)
						imin = 0;
					else {
						pxt = xs + tmax*(x_diff);
						imin = static_cast<int>(ceil((pxt - bx) / d));
					}
					tx0 = (bx + static_cast<double>(imax) * d - xs) / (x_diff);
					iu = -1;
				}


				if (ys < yd) {
					if (tmin == tymin)
						jmin = 1;
					else {
						pyt = ys + tmin*(y_diff);
						jmin = static_cast<int>(ceil((pyt - by) / d));
					}
					if (tmax == tymax)
						jmax = Ny;
					else {
						pyt = ys + tmax*(y_diff);
						jmax = static_cast<int>(floor((pyt - by) / d));
					}
					ty0 = (by + static_cast<double>(jmin) * d - ys) / (y_diff);
					ju = 1;

				}
				else if (ys > yd) {
					if (tmin == tymin)
						jmax = Ny - 1;
					else {
						pyt = ys + tmin*(y_diff);
						jmax = static_cast<int>(floor((pyt - by) / d));
					}
					if (tmax == tymax)
						jmin = 0;
					else {
						pyt = ys + tmax*(y_diff);
						jmin = static_cast<int>(ceil((pyt - by) / d));
					}
					ty0 = (by + static_cast<double>(jmax) * d - ys) / (y_diff);
					ju = -1;
				}

				int tempi, tempj, tempk;

				double pt = ((min(tx0, ty0) + tmin) / 2.);

				tempi = static_cast<int>(floor(((xs + pt * x_diff) - bx) / d));
				tempj = static_cast<int>(floor(((ys + pt * y_diff) - by) / d));

				double txu = d / fabs(x_diff);
				double tyu = d / fabs(y_diff);

				double L = sqrt(x_diff*x_diff + y_diff*y_diff);

				tempk = z_loop * Ny * Nx;

				double tc = tmin;

				double temp = 0.;

				double apu;

				unsigned int nnn = 0;
				int apu_tempi;
				bool apu_var = false;
				double tx0_apu, tc_apu;

				for (int ii = 0; ii < Np; ii++) {

					indices[ii] = tempj*Nx + tempi + tempk;

					if (tx0 < ty0) {
						elements[ii] = (tx0 - tc) * L;

						tempi += iu;
						tc = tx0;
						tx0 += txu;

						temp += elements[ii];


					}
					else {

						elements[ii] = (ty0 - tc) * L;

						tempj += ju;
						tc = ty0;
						ty0 += tyu;

						temp += elements[ii];
					}
				}

				temp = 1. / temp;

				if (attenuation_correction) {

					double jelppi = 0.;

					for (unsigned int iii = 0; iii < Np; iii++) {
						jelppi += (elements[iii]) * -atten[indices[iii]];
					}
					temp = exp(jelppi) * temp;

				}

				for (unsigned int ii = 0; ii < Np; ii++) {
					elements[ii] = elements[ii] * (temp);
				}
				pass = true;
			}
		}
		else {

			if (fabs(y_diff) < 1e-8) {
				if (yd <= maxyy && yd >= minyy) {

					double tx0 = (bxf - xs) / (x_diff);
					double tz0 = (bzf - zs) / (z_diff);
					double txback = (bxb - xs) / (x_diff);
					double tzback = (bzb - zs) / (z_diff);

					double txmin = min(tx0, txback);
					double txmax = max(tx0, txback);
					double tzmin = min(tz0, tzback);
					double tzmax = max(tz0, tzback);

					double tmin = max(txmin, tzmin);
					double tmax = min(txmax, tzmax);

					int imin, imax, kmin, kmax;
					double pxt, pzt;

					int iu, ku;

					if (xs < xd) {
						if (tmin == txmin)
							imin = 1;
						else {
							pxt = xs + tmin*(x_diff);
							imin = static_cast<int>(ceil((pxt - bx) / d));
						}
						if (tmax == txmax)
							imax = Nx;
						else {
							pxt = xs + tmax*(x_diff);
							imax = static_cast<int>(floor((pxt - bx) / d));
						}
						tx0 = (bx + static_cast<double>(imin) * d - xs) / (x_diff);
						iu = 1;
					}
					else if (xs > xd) {
						if (tmin == txmin)
							imax = Nx - 1;
						else {
							pxt = xs + tmin*(x_diff);
							imax = static_cast<int>(floor((pxt - bx) / d));
						}
						if (tmax == txmax)
							imin = 0;
						else {
							pxt = xs + tmax*(x_diff);
							imin = static_cast<int>(ceil((pxt - bx) / d));
						}
						tx0 = (bx + static_cast<double>(imax) * d - xs) / (x_diff);
						iu = -1;
					}

					if (zs < zd) {
						if (tmin == tzmin)
							kmin = 1;
						else {
							pzt = zs + tmin*(z_diff);
							kmin = static_cast<int>(ceil((pzt - bz) / dz));
						}
						if (tmax == tzmax)
							kmax = Nz;
						else {
							pzt = zs + tmax*(z_diff);
							kmax = static_cast<int>(floor((pzt - bz) / dz));
						}
						tz0 = (bz + static_cast<double>(kmin) * dz - zs) / (z_diff);
						ku = 1;
					}
					else if (zs > zd) {
						if (tmin == tzmin)
							kmax = Nz - 1;
						else {
							pzt = zs + tmin*(z_diff);
							kmax = static_cast<int>(floor((pzt - bz) / dz));
						}
						if (tmax == tzmax)
							kmin = 0;
						else {
							pzt = zs + tmax*(z_diff);
							kmin = static_cast<int>(ceil((pzt - bz) / dz));
						}
						tz0 = (bz + static_cast<double>(kmax) * dz - zs) / (z_diff);
						ku = -1;
					}

					double L = sqrt((x_diff*x_diff + z_diff*z_diff));

					int tempi, tempj, tempk;
					double apu2, apu1;

					double pt = ((min(tx0, tz0) + tmin) / 2.);

					tempi = static_cast<int>(floor(((xs + pt * x_diff) - bx) / d));
					tempk = static_cast<int>(floor(((zs + pt * z_diff) - bz) / dz));

					double txu = d / fabs(x_diff);
					double tzu = dz / fabs(z_diff);

					for (unsigned int ii = 0; ii < Np; ii++) {
						apu1 = fabs(yy_vec[ii] - yd);
						if (ii > 0 && apu1 < apu2 || ii == 0) {
							tempj = ii;
							apu2 = apu1;
						}
					}

					tempj = tempj * Nx;

					double tc = tmin;

					double temp = 0.;

					double apu;

					unsigned int nnn = 0;
					int apu_tempi;
					bool apu_var = false;
					double tx0_apu, tc_apu;

					for (int ii = 0; ii < Np; ii++) {

						indices[ii] = tempj + tempi + Nx * Ny * tempk;
						if (tx0 < tz0) {

							elements[ii] = (tx0 - tc) * L;


							tempi += iu;
							tc = tx0;
							tx0 += txu;

							temp += elements[ii];
						}
						else {

							elements[ii] = (tz0 - tc) * L;

							tempk += ku;
							tc = tz0;
							tz0 += tzu;

							temp += elements[ii];
						}

					}

					temp = 1. / temp;

					if (attenuation_correction) {

						double jelppi = 0.;

						for (unsigned int iii = 0; iii < Np; iii++) {
							jelppi += elements[iii] * -atten[indices[iii]];
						}
						temp = exp(jelppi) * temp;
					}

					for (unsigned int ii = 0; ii < Np; ii++) {
						elements[ii] = (elements[ii] * temp);
					}
					pass = true;

				}
			}
			else if (fabs(x_diff) < 1e-8) {
				if (xd <= maxxx && xd >= minxx) {

					double ty0 = (byf - ys) / (y_diff);
					double tyback = (byb - ys) / (y_diff);
					double tz0 = (bzf - zs) / (z_diff);
					double tzback = (bzb - zs) / (z_diff);

					double tzmin = min(tz0, tzback);
					double tzmax = max(tz0, tzback);
					double tymin = min(ty0, tyback);
					double tymax = max(ty0, tyback);

					double tmin = max(tymin, tzmin);
					double tmax = min(tymax, tzmax);

					int jmin, jmax, kmin, kmax;
					double pyt, pzt;

					int ku, ju;

					if (ys < yd) {
						if (tmin == tymin)
							jmin = 1;
						else {
							pyt = ys + tmin*(y_diff);
							jmin = static_cast<int>(ceil((pyt - by) / d));
						}
						if (tmax == tymax)
							jmax = Ny;
						else {
							pyt = ys + tmax*(y_diff);
							jmax = static_cast<int>(floor((pyt - by) / d));
						}
						ty0 = (by + static_cast<double>(jmin) * d - ys) / (y_diff);
						ju = 1;
					}
					else if (ys > yd) {
						if (tmin == tymin)
							jmax = Ny - 1;
						else {
							pyt = ys + tmin*(y_diff);
							jmax = static_cast<int>(floor((pyt - by) / d));
						}
						if (tmax == tymax)
							jmin = 0;
						else {
							pyt = ys + tmax*(y_diff);
							jmin = static_cast<int>(ceil((pyt - by) / d));
						}
						ty0 = (by + static_cast<double>(jmax) * d - ys) / (y_diff);
						ju = -1;
					}


					if (zs < zd) {
						if (tmin == tzmin)
							kmin = 1;
						else {
							pzt = zs + tmin*(z_diff);
							kmin = static_cast<int>(ceil((pzt - bz) / dz));
						}
						if (tmax == tzmax)
							kmax = Nz;
						else {
							pzt = zs + tmax*(z_diff);
							kmax = static_cast<int>(floor((pzt - bz) / dz));
						}
						tz0 = (bz + static_cast<double>(kmin) * dz - zs) / (z_diff);
						ku = 1;
					}
					else if (zs > zd) {
						if (tmin == tzmin)
							kmax = Nz - 1;
						else {
							pzt = zs + tmin*(z_diff);
							kmax = static_cast<int>(floor((pzt - bz) / dz));
						}
						if (tmax == tzmax)
							kmin = 0;
						else {
							pzt = zs + tmax*(z_diff);
							kmin = static_cast<int>(ceil((pzt - bz) / dz));
						}
						tz0 = (bz + static_cast<double>(kmax) * dz - zs) / (z_diff);
						ku = -1;
					}

					double L = sqrt((y_diff*y_diff + z_diff*z_diff));

					int tempi, tempj, tempk;
					double apu2, apu1;

					double pt = ((min(tz0, ty0) + tmin) / 2.);

					tempk = static_cast<int>(floor(((zs + pt * z_diff) - bz) / dz));
					tempj = static_cast<int>(floor(((ys + pt * y_diff) - by) / d));

					double tzu = dz / fabs(z_diff);
					double tyu = d / fabs(y_diff);

					double tc = tmin;

					double temp = 0.;

					double apu;

					for (int ii = 0; ii < Nx; ii++) {
						apu1 = fabs(xx_vec[ii] - xd);
						if (ii > 0 && apu1 < apu2 || ii == 0) {
							tempi = ii;
							apu2 = apu1;
						}
					}

					unsigned int nnn = 0;
					int apu_tempj;
					bool apu_var = false;
					double ty0_apu, tc_apu;

					for (int ii = 0; ii < Np; ii++) {
						indices[ii] = tempj * Nx + tempi + Nx * Ny * tempk;
						if (tz0 < ty0) {

							elements[ii] = (tz0 - tc) * L;

							tempk += ku;
							tc = tz0;
							tz0 += tzu;

							temp += elements[ii];

						}
						else {

							elements[ii] = (ty0 - tc) * L;

							tempj += ju;
							tc = ty0;
							ty0 += tyu;

							temp += elements[ii];
						}
					}

					temp = 1. / temp;

					if (attenuation_correction) {

						double jelppi = 0.;

						for (unsigned int iii = 0; iii < Np; iii++) {
							jelppi += elements[iii] * -atten[indices[iii]];
						}
						temp = exp(jelppi) * temp;
					}

					for (unsigned int ii = 0; ii < Np; ii++) {
						elements[ii] = (elements[ii] * temp);
					}
					pass = true;

				}
			}
			else {

				double tx0 = (bxf - xs) / (x_diff);
				double tz0 = (bzf - zs) / (z_diff);
				double txback = (bxb - xs) / (x_diff);
				double tzback = (bzb - zs) / (z_diff);
				double ty0 = (byf - ys) / (y_diff);
				double tyback = (byb - ys) / (y_diff);

				double txmin = min(tx0, txback);
				double txmax = max(tx0, txback);
				double tymin = min(ty0, tyback);
				double tymax = max(ty0, tyback);
				double tzmin = min(tz0, tzback);
				double tzmax = max(tz0, tzback);

				double tmin = max(max(txmin, tzmin), tymin);
				double tmax = min(min(txmax, tzmax), tymax);

				int imin, imax, jmin, jmax, kmin, kmax;
				double pxt, pyt, pzt;

				int iu, ju, ku;


				if (xs < xd) {
					if (tmin == txmin)
						imin = 1;
					else {
						pxt = xs + tmin*(x_diff);
						imin = static_cast<int>(ceil((pxt - bx) / d));
					}
					if (tmax == txmax)
						imax = Nx;
					else {
						pxt = xs + tmax*(x_diff);
						imax = static_cast<int>(floor((pxt - bx) / d));
					}
					tx0 = (bx + static_cast<double>(imin) * d - xs) / (x_diff);
					iu = 1;
				}
				else if (xs > xd) {
					if (tmin == txmin)
						imax = Nx - 1;
					else {
						pxt = xs + tmin*(x_diff);
						imax = static_cast<int>(floor((pxt - bx) / d));
					}
					if (tmax == txmax)
						imin = 0;
					else {
						pxt = xs + tmax*(x_diff);
						imin = static_cast<int>(ceil((pxt - bx) / d));
					}
					tx0 = (bx + static_cast<double>(imax) * d - xs) / (x_diff);
					iu = -1;
				}

				if (ys < yd) {
					if (tmin == tymin)
						jmin = 1;
					else {
						pyt = ys + tmin*(y_diff);
						jmin = static_cast<int>(ceil((pyt - by) / d));
					}
					if (tmax == tymax)
						jmax = Ny;
					else {
						pyt = ys + tmax*(y_diff);
						jmax = static_cast<int>(floor((pyt - by) / d));
					}
					ty0 = (by + static_cast<double>(jmin) * d - ys) / (y_diff);
					ju = 1;
				}
				else if (ys > yd) {
					if (tmin == tymin)
						jmax = Ny - 1;
					else {
						pyt = ys + tmin*(y_diff);
						jmax = static_cast<int>(floor((pyt - by) / d));
					}
					if (tmax == tymax)
						jmin = 0;
					else {
						pyt = ys + tmax*(y_diff);
						jmin = static_cast<int>(ceil((pyt - by) / d));
					}

					ty0 = (by + static_cast<double>(jmax) * d - ys) / (y_diff);
					ju = -1;
				}

				if (zs < zd) {
					if (tmin == tzmin)
						kmin = 1;
					else {
						pzt = zs + tmin*(z_diff);
						kmin = static_cast<int>(ceil((pzt - bz) / dz));
					}
					if (tmax == tzmax)
						kmax = Nz;
					else {
						pzt = zs + tmax*(z_diff);
						kmax = static_cast<int>(floor((pzt - bz) / dz));
					}

					tz0 = (bz + static_cast<double>(kmin) * dz - zs) / (z_diff);
					ku = 1;
				}
				else if (zs > zd) {
					if (tmin == tzmin)
						kmax = Nz - 1;
					else {
						pzt = zs + tmin*(z_diff);
						kmax = static_cast<int>(floor((pzt - bz) / dz));
					}
					if (tmax == tzmax)
						kmin = 0;
					else {
						pzt = zs + tmax*(z_diff);
						kmin = static_cast<int>(ceil((pzt - bz) / dz));
					}
					tz0 = (bz + static_cast<double>(kmax) * dz - zs) / (z_diff);
					ku = -1;
				}

				double L = sqrt(x_diff*x_diff + z_diff*z_diff + y_diff*y_diff);

				int tempi, tempj, tempk;

				double pt = ((min(min(tz0, ty0), tx0) + tmin) / 2.);

				tempk = static_cast<int>(floor(((zs + pt * z_diff) - bz) / dz));
				tempj = static_cast<int>(floor(((ys + pt * y_diff) - by) / d));
				tempi = static_cast<int>(floor(((xs + pt * x_diff) - bx) / d));

				double tzu = dz / fabs(z_diff);
				double tyu = d / fabs(y_diff);
				double txu = d / fabs(x_diff);

				double tc = tmin;

				double temp = 0.;

				double apu;

				unsigned int nnn = 0;
				int apu_tempi, apu_tempj;
				bool apu_varx = false, apu_vary = false, apu_varx2 = false;
				double tx0_apu, tc_apu, ty0_apu;
				bool final_x = false, final_y = false, first_x = false, first_y = false;
				unsigned int nnx = 0;
				unsigned int nny = 0;
				unsigned int nyy = 0;
				unsigned int nxx = 0;

				for (int ii = 0; ii < Np; ii++) {


					indices[ii] = tempj * Nx + tempi + Nx * Ny * tempk;
					if (tz0 < ty0 && tz0 < tx0) {

						elements[ii] = (tz0 - tc) * L;

						tempk += ku;
						tc = tz0;
						tz0 += tzu;

						temp += elements[ii];
					}
					else if (ty0 < tx0) {
						elements[ii] = (ty0 - tc) * L;

						tempj += ju;
						tc = ty0;
						ty0 += tyu;

						temp += elements[ii];
					}
					else {
						elements[ii] = (tx0 - tc) * L;

						tempi += iu;
						tc = tx0;
						tx0 += txu;

						temp += elements[ii];
					}

				}

				temp = 1. / temp;

				if (attenuation_correction) {

					double jelppi = 0.;

					for (unsigned int iii = 0; iii < Np; iii++) {
						jelppi += elements[iii] * -atten[indices[iii]];
					}
					temp = exp(jelppi) * temp;

				}

				for (unsigned int ii = 0; ii < Np; ii++) {
					elements[ii] = (elements[ii] * temp);
				}
				pass = true;
			}
		}
		double ax = 0.;
		double yax = 0.;
		if (pass) {
			if (Sino[oo] != 0.) {
				for (int k = 0; k < Np; k++) {
					ax += (elements[k] * osem_apu[indices[k]]);
				}
				ax += epps;
				yax = Sino[oo] / ax;
				for (int k = 0; k < Np; k++) {
					rhs[indices[k]] += (elements[k] * yax);
					Summ[indices[k]] += elements[k];
				}
			}
			else {
				for (int k = 0; k < Np; k++) {
					Summ[indices[k]] += elements[k];
				}
			}
		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray*prhs[])

{

	if (nrhs != 32)
		mexErrMsgTxt("Invalid number of input arguments.  There must be exactly 32.");

	if (nlhs != 2)
		mexErrMsgTxt("Invalid number of output arguments.  There must be exactly two.");


	int Ny = (int)mxGetScalar(prhs[0]);

	int Nx = (int)mxGetScalar(prhs[1]);

	int Nz = (int)mxGetScalar(prhs[2]);

	double d = (double)mxGetScalar(prhs[3]);

	double dz = (double)mxGetScalar(prhs[4]);

	double by = (double)mxGetScalar(prhs[5]);

	double bx = (double)mxGetScalar(prhs[6]);

	double bz = (double)mxGetScalar(prhs[7]);

	double *z_det = (double*)mxGetData(prhs[8]);

	double *x = (double*)mxGetData(prhs[9]);

	double *y = (double*)mxGetData(prhs[10]);

	double *iij = (double*)mxGetData(prhs[11]);

	double *jji = (double*)mxGetData(prhs[12]);

	double *kkj = (double*)mxGetData(prhs[13]);

	double *yy = (double*)mxGetData(prhs[14]);

	double *xx = (double*)mxGetData(prhs[15]);

	int NSinos = (int)mxGetScalar(prhs[16]);

	int NSlices = (int)mxGetScalar(prhs[17]);

	vector<double> iij_vec(iij, iij + mxGetNumberOfElements(prhs[11]));

	vector<double> jjk_vec(jji, jji + mxGetNumberOfElements(prhs[12]));

	vector<double> kkj_vec(kkj, kkj + mxGetNumberOfElements(prhs[13]));

	vector<double> yy_vec(yy, yy + mxGetNumberOfElements(prhs[14]) - 1);

	vector<double> xx_vec(xx, xx + mxGetNumberOfElements(prhs[15]) - 1);

	int size_x = (int)mxGetScalar(prhs[18]);

	double zmax = (double)mxGetScalar(prhs[19]);

	double *atten = (double*)mxGetData(prhs[20]);

	int pituus = (int)mxGetScalar(prhs[21]);

	bool attenuation_correction = (bool)mxGetScalar(prhs[22]);

	uint16_t* lor1 = (uint16_t*)mxGetData(prhs[23]);

	uint32_t* xy_index = (uint32_t*)mxGetData(prhs[24]);

	uint32_t* z_index = (uint32_t*)mxGetData(prhs[25]);

	int TotSinos = (int)mxGetScalar(prhs[26]);

	int maksimi = (int)mxGetScalar(prhs[27]);

	double epps = (double)mxGetScalar(prhs[28]);

	double* Sino = (double*)mxGetData(prhs[29]);

	double* osem_apu = (double*)mxGetData(prhs[30]);

	bool verbose = (bool)mxGetScalar(prhs[31]);

	size_t loop_var_par = pituus;

	size_t N = Nx * Ny * Nz;

	double maxyy = *max_element(yy, yy + Ny + 1);
	double minyy = *min_element(yy, yy + Ny + 1);

	double maxxx = *max_element(xx, xx + Nx + 1);
	double minxx = *min_element(xx, xx + Nx + 1);

	plhs[0] = mxCreateNumericMatrix(N, 1, mxDOUBLE_CLASS, mxREAL);

	double* Summ = (double*)mxGetData(plhs[0]);

	plhs[1] = mxCreateNumericMatrix(N, 1, mxDOUBLE_CLASS, mxREAL);

	double* rhs = (double*)mxGetData(plhs[1]);

	clock_t time = clock();

	improved_siddon(loop_var_par, size_x, zmax, Summ, rhs, maxyy, minyy, maxxx, minxx, xx_vec, iij_vec, jjk_vec, kkj_vec, yy_vec, atten, x, y, z_det, NSlices, Nx, Ny, Nz, d, dz,
		bx, by, bz, attenuation_correction, lor1, xy_index, z_index, TotSinos, maksimi, epps, Sino, osem_apu, N);

	time = clock() - time;

	if (verbose) {
		mexPrintf("Improved Siddon took %f seconds\n", ((float)time) / CLOCKS_PER_SEC);
		mexEvalString("pause(.001);");
	}

	return;
}