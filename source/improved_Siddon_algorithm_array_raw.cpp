#include "mex.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <time.h>
#include <cstdint>


using namespace std;

void improved_siddon(const size_t loop_var_par, const int size_x, const double zmax, mwIndex* indices, double* elements, mwIndex* lor, const double maxyy, 
	const double minyy,	const double maxxx, const double minxx, const vector<double>& xx_vec, const vector<double>& iij_vec, const vector<double>& jjk_vec, 
	const vector<double>& kkj_vec, const vector<double>& yy_vec, const double* atten, const double* x, const double* y, const double* z_det, const int Nx, 
	const int Ny, const int Nz, const double d, const double dz, const double bx, const double by, const double bz, const bool attenuation_correction, 
	const uint16_t* lor1, const unsigned int* lor2, const uint16_t *L, const int *pseudos, const int pRows, const int det_per_ring) {

	int oo = -1;
	int ll = 0;
	int ps;

	const double bxf = bx + iij_vec.front() * d;
	const double byf = by + jjk_vec.front() * d;
	const double bxb = bx + iij_vec.back() * d;
	const double byb = by + jjk_vec.back() * d;
	const double bzf = bz + kkj_vec.front() * dz;
	const double bzb = bz + kkj_vec.back() * dz;

	lor[0] = 0;

	for (int lo = 1; lo <= loop_var_par; lo++) {
			
		oo++;

		int detektorit1 = static_cast<int>(L[ll++]) - 1;
		int detektorit2 = static_cast<int>(L[ll++]) - 1;

		int loop1 = ((detektorit1) / det_per_ring);
		int loop2 = ((detektorit2) / det_per_ring);

		double zs, zd;

		if (loop1 == loop2) {
			zs = z_det[loop1];
			zd = zs;
		}
		else {
			zs = z_det[loop1];
			zd = z_det[loop2];
		}

		

		if (loop1 >= pseudos[0]) {
			ps = 1;
			for (int kk = 0; kk < pRows; kk++) {
				if (kk + 1 < pRows) {
					if (loop1 >= pseudos[kk] && loop1 < pseudos[kk + 1]) {
						zs = z_det[loop1 + ps];
						break;
					}
					else
						ps++;
				}
				else {
					if (loop1 >= pseudos[kk])
						zs = z_det[loop1 + ps];
				}
			}
		}

		if (loop2 >= pseudos[0]) {
			ps = 1;

			for (int kk = 0; kk < pRows; kk++) {
				if (kk + 1 < pRows) {
					if (loop2 >= pseudos[kk] && loop2 < pseudos[kk + 1]) {
						zd = z_det[loop2 + ps];
						break;
					}
					else
						ps++;
				}
				else {
					if (loop2 >= pseudos[kk])
						zd = z_det[loop2 + ps];
				}
			}
		}

		double xs = x[detektorit1 - det_per_ring*(loop1)];
		double xd = x[detektorit2 - det_per_ring*(loop2)];
		double ys = y[detektorit1 - det_per_ring*(loop1)];
		double yd = y[detektorit2 - det_per_ring*(loop2)];
		

		double y_diff = (yd - ys);
		double x_diff = (xd - xs);
		double z_diff = (zd - zs);

		unsigned int Np = static_cast<unsigned int>(lor1[oo]);
		unsigned int N2 = lor2[oo];
		

		if (fabs(z_diff) < 1e-8) {

			int z_loop = static_cast<int>((zs / zmax)*(Nz - 1));

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
						indices[N2 + ii] = tempk + ii;
						elements[N2 + ii] = fabs(element);
					}
					lor[lo] = lor2[lo];
					continue;
				}
				continue;
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
						indices[N2 + ii] = tempk + ii * Ny;
						elements[N2 + ii] = fabs(element);
					}
					lor[lo] = lor2[lo];
					continue;
				}

				continue;
			}
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

			for (unsigned int ii = 0; ii < Np; ii++) {

				if (tx0 < ty0) {

					if (iu > 0) {

						apu = (tx0 - tc) * L;
						elements[N2 + ii] = (apu);
						indices[N2 + ii] = tempj*Nx + tempi + tempk;

						tempi += iu;
						tc = tx0;
						tx0 += txu;

						temp += apu;
					}
					else {
						if (nnn == 0) {
							apu_var = true;
						}
						tempi += iu;
						tx0 += txu;
						nnn++;
					}


				}
				else {

					if (apu_var) {
							apu_tempi = tempi;
							tx0_apu = tx0;
							tc_apu = tx0 - txu;
						for (int nl = static_cast<int>(nnn); nl >= 0; nl--) {
							indices[N2 + ii - nl] = tempj*Nx + apu_tempi + tempk;
							if (nl == static_cast<int>(nnn)) {
								apu = (ty0 - tc_apu) * L;
							}
							else if (nl == 0)
								apu = (tx0_apu - tc) * L;
							else {
								apu = (tx0_apu - tc_apu) * L;
							}
							elements[N2 + ii - nl] = (apu);
							temp += apu;

							apu_tempi -= iu;
							tx0_apu -= txu;
							tc_apu -= txu;
						}
						apu_var = false;
						nnn = 0;

						tempj += ju;
						tc = ty0;
						ty0 += tyu;
					}
					else {
						apu = (ty0 - tc) * L;
						elements[N2 + ii] = (apu);
						indices[N2 + ii] = tempj*Nx + tempi + tempk;

						tempj += ju;
						tc = ty0;
						ty0 += tyu;

						temp += apu;
					}
				}

				if (ii == Np - 1 && apu_var) {
						apu_tempi = tempi - iu;
						tx0_apu = tx0 - txu;
						tc_apu = tx0_apu - txu;
						for (int nl = static_cast<int>(nnn) - 1; nl >= 0; nl--) {
							indices[N2 + ii - nl] = tempj*Nx + apu_tempi + tempk;
							if (nl == 0)
								apu = (tx0_apu - tc) * L;
							else {
								apu = (tx0_apu - tc_apu) * L;
							}
							elements[N2 + ii - nl] = (apu);
							temp += apu;
							if (nnn > 1) {
								apu_tempi -= iu;
								tx0_apu -= txu;
								tc_apu -= txu;
							}
						}
					apu_var = false;
					nnn = 0;
				}

			}

			temp = 1. / temp;

			if (attenuation_correction) {

				double jelppi = 0.;

				for (unsigned int iii = 0; iii < Np; iii++) {
					jelppi += (elements[N2 + iii]) * -atten[indices[N2 + iii]];
				}
				temp = exp(jelppi) * temp;

			}


			for (unsigned int ii = 0; ii < Np; ii++) {
				elements[N2 + ii] = fabs(elements[N2 + ii]) * (temp);
			}
			lor[lo] = lor2[lo];
			continue;
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

					for (unsigned int ii = 0; ii < Np; ii++) {

						if (tx0 < tz0) {

							if (iu > 0) {

								apu = (tx0 - tc) * L;
								elements[N2 + ii] = (apu);
								indices[N2 + ii] = tempj + tempi + Nx * Ny * tempk;

								tempi += iu;
								tc = tx0;
								tx0 += txu;

								temp += apu;
							}
							else {
								if (nnn == 0) {
									apu_var = true;
								}
								tempi += iu;
								tx0 += txu;
								nnn++;
							}


						}
						else {

							if (apu_var) {
								apu_tempi = tempi;
								tx0_apu = tx0;
								tc_apu = tx0 - txu;
								for (int nl = static_cast<int>(nnn); nl >= 0; nl--) {
									indices[N2 + ii - nl] = tempj + apu_tempi + Nx * Ny * tempk;
									if (nl == static_cast<int>(nnn)) {
										apu = (tz0 - tc_apu) * L;
									}
									else if (nl == 0)
										apu = (tx0_apu - tc) * L;
									else {
										apu = (tx0_apu - tc_apu) * L;
									}
									elements[N2 + ii - nl] = (apu);
									temp += apu;

									apu_tempi -= iu;
									tx0_apu -= txu;
									tc_apu -= txu;
								}
								apu_var = false;
								nnn = 0;

								tempk += ku;
								tc = tz0;
								tz0 += tzu;
							}
							else {
								apu = (tz0 - tc) * L;
								elements[N2 + ii] = (apu);
								indices[N2 + ii] = tempj + tempi + Nx * Ny * tempk;

								tempk += ku;
								tc = tz0;
								tz0 += tzu;

								temp += apu;
							}
						}

						if (ii == Np - 1 && apu_var) {
							apu_tempi = tempi - iu;
							tx0_apu = tx0 - txu;
							tc_apu = tx0_apu - txu;
							for (int nl = static_cast<int>(nnn) - 1; nl >= 0; nl--) {
								indices[N2 + ii - nl] = tempj + apu_tempi + Nx * Ny * tempk;
								if (nl == 0)
									apu = (tx0_apu - tc) * L;
								else {
									apu = (tx0_apu - tc_apu) * L;
								}
								elements[N2 + ii - nl] = (apu);
								temp += apu;
								if (nnn > 1) {
									apu_tempi -= iu;
									tx0_apu -= txu;
									tc_apu -= txu;
								}
							}
							apu_var = false;
							nnn = 0;
						}
					}

					temp = 1. / temp;

					if (attenuation_correction) {

						double jelppi = 0.;

						for (unsigned int iii = 0; iii < Np; iii++) {
							jelppi += elements[N2 + iii] * -atten[indices[N2 + iii]];
						}
						temp = exp(jelppi) * temp;
					}

					for (unsigned int ii = 0; ii < Np; ii++) {
						elements[N2 + ii] = (fabs(elements[N2 + ii]) * temp);
					}
					lor[lo] = lor2[lo];
					continue;

				}
				continue;
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

					for (unsigned int ii = 0; ii < Np; ii++) {

						if (tz0 < ty0) {
							if (apu_var) {
								apu_tempj = tempj;
								ty0_apu = ty0;
								tc_apu = ty0 - tyu;
								for (int nl = static_cast<int>(nnn); nl >= 0; nl--) {
									indices[N2 + ii - nl] = apu_tempj * Nx + tempi + Nx * Ny * tempk;
									if (nl == static_cast<int>(nnn)) {
										apu = (tz0 - tc_apu) * L;
									}
									else if (nl == 0)
										apu = (ty0_apu - tc) * L;
									else {
										apu = (ty0_apu - tc_apu) * L;
									}
									elements[N2 + ii - nl] = (apu);
									temp += apu;

									apu_tempj -= ju;
									ty0_apu -= tyu;
									tc_apu -= tyu;
								}
								apu_var = false;
								nnn = 0;

								tempk += ku;
								tc = tz0;
								tz0 += tzu;
							}
							else {
								apu = (tz0 - tc) * L;
								elements[N2 + ii] = (apu);
								indices[N2 + ii] = tempj * Nx + tempi + Nx * Ny * tempk;

								tempk += ku;
								tc = tz0;
								tz0 += tzu;

								temp += apu;
							}

							

						}
						else {

							if (ju > 0) {
								apu = (ty0 - tc) * L;
								elements[N2 + ii] = (apu);
								indices[N2 + ii] = tempj * Nx + tempi + Nx * Ny * tempk;

								tempj += ju;
								tc = ty0;
								ty0 += tyu;

								temp += apu;
							}
							else {
								if (nnn == 0) {
									apu_var = true;
								}
								tempj += ju;
								//tc = ty0;
								ty0 += tyu;
								nnn++;
							}
						}

						if (ii == Np - 1 && apu_var) {
							apu_tempj = tempj - ju;
							ty0_apu = ty0 - tyu;
							tc_apu = ty0_apu - tyu;
							for (int nl = static_cast<int>(nnn) - 1; nl >= 0; nl--) {
								indices[N2 + ii - nl] = apu_tempj * Nx + tempi + Nx * Ny * tempk;
								if (nl == 0)
									apu = (ty0_apu - tc) * L;
								else {
									apu = (ty0_apu - tc_apu) * L;
								}
								elements[N2 + ii - nl] = (apu);
								temp += apu;
								if (nnn > 1) {
									apu_tempj -= ju;
									ty0_apu -= tyu;
									tc_apu -= tyu;
								}
							}
							apu_var = false;
							nnn = 0;
						}
					}

					temp = 1. / temp;

					if (attenuation_correction) {

						double jelppi = 0.;

						for (unsigned int iii = 0; iii < Np; iii++) {
							jelppi += elements[N2 + iii] * -atten[indices[N2 + iii]];
						}
						temp = exp(jelppi) * temp;
					}

					for (unsigned int ii = 0; ii < Np; ii++) {
						elements[N2 + ii] = (fabs(elements[N2 + ii]) * temp);
					}
					lor[lo] = lor2[lo];
					continue;

				}
				continue;
			}

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

			for (unsigned int ii = 0; ii < Np; ii++) {


				if (tz0 < ty0 && tz0 < tx0) {

					if (apu_varx && apu_vary == false) {
						apu_tempi = tempi;
						tx0_apu = tx0;
						tc_apu = tx0 - txu;
						for (int nl = static_cast<int>(nnn); nl >= 0; nl--) {
							indices[N2 + ii - nl] = tempj*Nx + apu_tempi + Nx * Ny * tempk;
							if (nl == static_cast<int>(nnn))
								apu = (tz0 - tc_apu) * L;
							else if (nl == 0)
								apu = (tx0_apu - tc) * L;
							else {
								apu = (tx0_apu - tc_apu) * L;
							}
							elements[N2 + ii - nl] = (apu);
							temp += apu;

							apu_tempi -= iu;
							tx0_apu -= txu;
							tc_apu -= txu;
						}
						apu_varx = false;
						nnn = 0;
					}
					else if (apu_vary && apu_varx == false) {
						if (apu_varx2) {
							if (final_x) {
								apu_tempi = tempi;
								tx0_apu = tx0;
								apu_tempj = tempj;
								ty0_apu = ty0;
							}
							else {
								apu_tempj = tempj - ju;
								ty0_apu = ty0 - tyu;
								apu_tempi = tempi;
								tx0_apu = tx0;
								tc_apu = ty0_apu;
								apu = (tz0 - tc_apu) * L;
								elements[N2 + ii - nnn - nnx] = (apu);
								indices[N2 + ii - nnn - nnx] = tempj*Nx + tempi + Nx * Ny * tempk;
								temp += apu;
								nnn--;
							}
							
							double tx0_apu2 = tx0_apu;
							double tc_apu2 = tc;
							int apu_tempi2 = tempi;

							for (int nl = static_cast<int>(nnn); nl >= 0; nl--) {
								unsigned int nnnx = 0;
								while (tx0_apu2 - txu >= ty0_apu - tyu && nnnx < nnx) {
									nnnx++;
									tx0_apu2 -= txu;
									apu_tempi2 -= iu;
								}
								if (nl == 0 && first_x)
									tc_apu2 = tc;
								else if (nnnx > 0)
									tc_apu2 = ty0_apu - tyu;

								double tx0_apu3 = tx0_apu2;
								int apu_tempi3 = apu_tempi2;
								if (nnnx > 0) {
									for (unsigned int nlx = 0; nlx < (nnnx); nlx++) {
										indices[N2 + ii - nl - (nnx - nlx)] = apu_tempj*Nx + apu_tempi2 + Nx * Ny * tempk;
										apu = (tx0_apu2 - tc_apu2) * L;
										elements[N2 + ii - nl - (nnx - nlx)] = (apu);
										temp += apu;
										tc_apu2 = tx0_apu2;
										if (nlx < (nnnx) - 1) {
											apu_tempi2 += iu;
											tx0_apu2 += txu;
										}
									}
									nnx -= nnnx;
									tx0_apu2 = tx0_apu3;
									apu_tempi2 = apu_tempi3;
									apu_tempi3 += nnnx * iu;
								}
								if (nnnx > 0)
									indices[N2 + ii - nl - nnx] = apu_tempj*Nx + apu_tempi3 + Nx * Ny * tempk;
								else
									indices[N2 + ii - nl - nnx] = apu_tempj*Nx + apu_tempi2 + Nx * Ny * tempk;
								if (nl == static_cast<int>(nnn)) {
									if (nnnx)
										tc_apu2 = tx0_apu - txu;
									else
										tc_apu2 = ty0_apu - tyu;
								}
								else if (nl > 0 && nnnx > 0 || nl == 0 && nnnx > 0) {
								}
								else if (nl > 0) {
									tc_apu2 = ty0_apu - tyu;
								}
								else if (first_y) {
									tc_apu2 = tc;
								}
								if (nl == static_cast<int>(nnn) && final_x)
									apu = (tz0 - tc_apu2) * L;
								else
									apu = (ty0_apu - tc_apu2) * L;
								elements[N2 + ii - nl - nnx] = (apu);
								temp += apu;

								apu_tempj -= ju;
								ty0_apu -= tyu;
							}
							apu_vary = false;
							nnn = 0;
							apu_varx2 = false;
							nnx = 0;
							first_y = false;
							first_x = false;
						}
						else {
							apu_tempj = tempj;
							ty0_apu = ty0;
							tc_apu = ty0 - tyu;
							for (int nl = static_cast<int>(nnn); nl >= 0; nl--) {
								indices[N2 + ii - nl] = apu_tempj*Nx + tempi + Nx * Ny * tempk;
								if (nl == static_cast<int>(nnn)) {
									apu = (tz0 - tc_apu) * L;
								}
								else if (nl == 0)
									apu = (ty0_apu - tc) * L;
								else {
									apu = (ty0_apu - tc_apu) * L;
								}
								elements[N2 + ii - nl] = (apu);
								temp += apu;

								apu_tempj -= ju;
								ty0_apu -= tyu;
								tc_apu -= tyu;
							}
							apu_vary = false;
							nnn = 0;
							first_y = false;
						}
					}
					else if (apu_varx && apu_vary) {
						bool apu_bool = false;
						if (final_x) {
							apu_tempi = tempi;
							tx0_apu = tx0;
							apu_tempj = tempj;
							ty0_apu = ty0;
							tc_apu = tx0_apu - txu;
						}
						else {
							apu_tempj = tempj;
							ty0_apu = ty0;
							apu_tempi = tempi;
							tx0_apu = tx0;
							tc_apu = ty0_apu - tyu;
						}
						apu = (tz0 - tc_apu) * L;
						elements[N2 + ii - nnn] = (apu);
						indices[N2 + ii - nnn] = tempj*Nx + tempi + Nx * Ny * tempk;
						temp += apu;
						nnn--;
						for (int nl = static_cast<int>(nnn); nl >= 0; nl--) {

							
							if (apu_bool) {
								tx0_apu -= txu;
								apu_tempi -= iu;
								if (tx0_apu - txu >= ty0_apu - tyu && nl > 0) {
									tc_apu = tx0_apu - txu;
								}
								else if (nl > 0) {
									tc_apu = ty0_apu - tyu;
									apu_bool = false;
								}
								else
									tc_apu = tc;
								apu = (tx0_apu - tc_apu) * L;
							}
							else {
								ty0_apu -= tyu;
								apu_tempj -= ju;
								if (tx0_apu - txu >= ty0_apu - tyu && nl > 0) {
									tc_apu = tx0_apu - txu;
									apu_bool = true;
								}
								else if (nl > 0) {
									tc_apu = ty0_apu - tyu;
								}
								else
									tc_apu = tc;
								apu = (ty0_apu - tc_apu) * L;
							}
							indices[N2 + ii - nl] = apu_tempj*Nx + apu_tempi + Nx * Ny * tempk;
							elements[N2 + ii - nl] = (apu);
							temp += apu;
						}
						apu_varx = false;
						apu_vary = false;
						nnn = 0;
						first_y = false;
					}
					else {
						double tc_apu2 = tc;
						if (apu_varx2) {
							apu_tempi = tempi;
							tx0_apu = tx0;
							for (unsigned int kk = 0; kk < (nnx); kk++) {
								tx0_apu -= txu;
								apu_tempi -= iu;
							}
							for (unsigned int nlx = 0; nlx < (nnx); nlx++) {
								indices[N2 + ii - (nnx - nlx)] = tempj*Nx + apu_tempi + Nx * Ny * tempk;
								apu_tempi += iu;
								apu = (tx0_apu - tc_apu2) * L;
								elements[N2 + ii - (nnx - nlx)] = (apu);
								temp += apu;
								tc_apu2 = tx0_apu;
								tx0_apu += txu;
							}
							apu_varx2 = false;
							nnx = 0;
							first_y = false;
							first_x = false;
						}
						indices[N2 + ii] = tempj*Nx + tempi + Nx * Ny * tempk;
						apu = ((tz0 - tc_apu2) * L);
						elements[N2 + ii] = apu;
						temp += apu;
					}
					tempk += ku;
					tc = tz0;
					tz0 += tzu;
				}
				else if (ty0 < tx0) {
					
					if (ju > 0) {

						if (apu_varx) {
							apu_tempi = tempi;
							tx0_apu = tx0;
							tc_apu = tx0 - txu;
							for (int nl = static_cast<int>(nnn); nl >= 0; nl--) {
								indices[N2 + ii - nl] = tempj*Nx + apu_tempi + Nx * Ny * tempk;
								if (nl == nnn)
									apu = (ty0 - tc_apu) * L;
								else if (nl == 0)
									apu = (tx0_apu - tc) * L;
								else {
									apu = (tx0_apu - tc_apu) * L;
								}
								elements[N2 + ii - nl] = (apu);
								temp += apu;

								apu_tempi -= iu;
								tx0_apu -= txu;
								tc_apu -= txu;
							}
							apu_varx = false;
							nnn = 0;
							tc = ty0;
						}
						else {

							indices[N2 + ii] = tempj*Nx + tempi + Nx * Ny * tempk;
							apu = ((ty0 - tc) * L);
							elements[N2 + ii] = apu;
							temp += apu;
							tc = ty0;
						}
					}
					else {
						if (nnn == 0 || apu_varx && apu_vary == false) {
							apu_vary = true;
						}
						final_y = true;
						final_x = false;
						nnn++;
						nyy++;
						if (first_y == false && first_x == false)
							first_y = true;
					}

					tempj += ju;
					ty0 += tyu;
				}
				else {

					if (iu > 0) {

						if (ju < 0 && ii < Np - 1) {
							apu_varx2 = true;
							nnx++;
							final_x = true;
							final_y = false;
							if (first_y == false && first_x == false)
								first_x = true;
						}
						else if (ii == Np - 1 && apu_varx2) {
							if (final_x) {
								tc_apu = tc;
								apu_tempi = tempi;
								tx0_apu = tx0;
								apu_tempj = tempj;
								ty0_apu = ty0;
							}
							else {
								apu_tempj = tempj - ju;
								ty0_apu = ty0 - tyu;
								apu_tempi = tempi;
								tx0_apu = tx0;
								apu = (tx0 - ty0_apu) * L;
								elements[N2 + ii - nnn - nnx] = (apu);
								indices[N2 + ii - nnn - nnx] = tempj*Nx + tempi + Nx * Ny * tempk;
								nnn--;
								temp += apu;
							}

							double tx0_apu2 = tx0_apu;
							double tc_apu2 = tc_apu;
							int apu_tempi2 = tempi;

							for (int nl = static_cast<int>(nnn); nl >= 0; nl--) {
								unsigned int nnnx = 0;
								bool nxn = false;
								while (tx0_apu2 - txu > ty0_apu - tyu && nnnx < nnx) {
									nnnx++;
									tx0_apu2 -= txu;
									apu_tempi2 -= iu;
								}
								if (nl == 0 && first_x) {
									tc_apu2 = tc;
								}
								else if (nnnx > 0)
									tc_apu2 = ty0_apu - tyu;
								double tx0_apu3 = tx0_apu2;
								int apu_tempi3 = apu_tempi2;
								if (nnnx > 0) {
									for (int nlx = 0; nlx < static_cast<int>(nnnx); nlx++) {
										indices[N2 + ii - nl - (nnx - nlx)] = apu_tempj*Nx + apu_tempi2 + Nx * Ny * tempk;
										apu = (tx0_apu2 - tc_apu2) * L;
										elements[N2 + ii - nl - (nnx - nlx)] = (apu);
										temp += apu;
										tc_apu2 = tx0_apu2;
										if (nlx < static_cast<int>(nnnx) -1) {
											apu_tempi2 += iu;
											tx0_apu2 += txu;
										}
									}
									nnx -= nnnx;
									nxn = true;
									tx0_apu2 = tx0_apu3;
									apu_tempi2 = apu_tempi3;
									apu_tempi3 += nnnx * iu;
								}
								if (nxn)
									indices[N2 + ii - nl - nnx] = apu_tempj*Nx + apu_tempi3 + Nx * Ny * tempk;
								else
									indices[N2 + ii - nl - nnx] = apu_tempj*Nx + apu_tempi2 + Nx * Ny * tempk;
								if (nnnx > 0 && nl > 0 || nl == 0 && nnnx > 0) {
								}
								else if (nl > 0) {
									tc_apu2 = ty0_apu - tyu;
								}
								else if (first_y) {
									tc_apu2 = tc;
								}
								if (nl == static_cast<int>(nnn) && final_x)
									apu = (tx0_apu - tc_apu2) * L;
								else
									apu = (ty0_apu - tc_apu2) * L;
								elements[N2 + ii - nl - nnx] = (apu);
								temp += apu;
								apu_tempj -= ju;
								ty0_apu -= tyu;
							}
							apu_vary = false;
							nnn = 0;
							apu_varx2 = false;
							nnx = 0;
							first_x = false;
							first_y = false;
						}
						else if (ii == Np - 1 && apu_varx2 == false && apu_vary == false) {
							indices[N2 + ii] = tempj*Nx + tempi + Nx * Ny * tempk;
							apu = ((tx0 - tc) * L);
							elements[N2 + ii] = apu;
							temp += apu;
						}
						else if (apu_vary && apu_varx2 == false) {
							apu_tempj = tempj;
							ty0_apu = ty0;
							tc_apu = ty0_apu - tyu;
							for (int nl = static_cast<int>(nnn); nl >= 0; nl--) {
								indices[N2 + ii - nl] = apu_tempj*Nx + tempi + Nx * Ny * tempk;
								if (nl == nnn) {
									apu = (tx0 - tc_apu) * L;
								}
								else if (nl == 0)
									apu = (ty0_apu - tc) * L;
								else {
									apu = (ty0_apu - tc_apu) * L;
								}
								elements[N2 + ii - nl] = (apu);
								temp += apu;

								apu_tempj -= ju;
								ty0_apu -= tyu;
								tc_apu -= tyu;
							}
							apu_vary = false;
							nnn = 0;
							tc = tx0;
							first_y = false;
						}
						else {
							indices[N2 + ii] = tempj*Nx + tempi + Nx * Ny * tempk;
							apu = ((tx0 - tc) * L);
							elements[N2 + ii] = apu;
							temp += apu;
							tc = tx0;
						}
					}
					else {
						if (nnn == 0 || apu_vary && apu_varx == false) {
							apu_varx = true;
						}
						final_y = false;
						final_x = true;
						nnn++;
						nxx++;
					}

					tempi += iu;
					tx0 += txu;
					
				}

				if (ii == Np - 1 && apu_varx && apu_vary == false) {
					apu_tempi = tempi - iu;
					tx0_apu = tx0 - txu;
					tc_apu = tx0_apu - txu;
					for (int nl = static_cast<int>(nnn) - 1; nl >= 0; nl--) {
						indices[N2 + ii - nl] = tempj*Nx + apu_tempi + Nx * Ny * tempk;
						if (nl == 0) {
							apu = (tx0_apu - tc) * L;
						}
						else {
							apu = (tx0_apu - tc_apu) * L;
						}
						elements[N2 + ii - nl] = (apu);
						temp += apu;

						if (nnn > 1) {
							apu_tempi -= iu;
							tx0_apu -= txu;
							tc_apu -= txu;
						}
					}
					apu_varx = false;
					nnn = 0;
				}
				else if (ii == Np - 1 && apu_vary && apu_varx == false) {
					if (apu_varx2) {
						apu_tempj = tempj - ju;
						ty0_apu = ty0 - tyu;
						apu_tempi = tempi;
						tx0_apu = tx0;

						double tx0_apu2 = tx0_apu;
						double tc_apu2 = tc_apu;
						int apu_tempi2 = tempi;

						for (int nl = static_cast<int>(nnn) - 1; nl >= 0; nl--) {
							unsigned int nnnx = 0;
							bool nxn = false;
							while (tx0_apu2 - txu > ty0_apu - tyu && nnnx < nnx) {
								nnnx++;
								tx0_apu2 -= txu;
								apu_tempi2 -= iu;
							}
							if (nl == 0 && first_x)
								tc_apu2 = tc;
							else if (nnnx > 0)
								tc_apu2 = ty0_apu - tyu;
							double tx0_apu3 = tx0_apu2;
							int apu_tempi3 = apu_tempi2;
							if (nnnx > 0) {
								for (int nlx = 0; nlx < static_cast<int>(nnnx); nlx++) {
									indices[N2 + ii - nl - (nnx - nlx)] = apu_tempj*Nx + apu_tempi2 + Nx * Ny * tempk;
									apu = (tx0_apu2 - tc_apu2) * L;
									elements[N2 + ii - nl - (nnx - nlx)] = (apu);
									temp += apu;
									tc_apu2 = tx0_apu2;
									if (nlx < static_cast<int>(nnnx) -1) {
										apu_tempi2 += iu;
										tx0_apu2 += txu;
									}
								}
								nnx -= nnnx;
								nxn = true;
								tx0_apu2 = tx0_apu3;
								apu_tempi2 = apu_tempi3;
								apu_tempi3 += nnnx * iu;
							}
							if (nxn) {
								indices[N2 + ii - nl - nnx] = apu_tempj*Nx + apu_tempi3 + Nx * Ny * tempk;
							}
							else {
								indices[N2 + ii - nl - nnx] = apu_tempj*Nx + apu_tempi2 + Nx * Ny * tempk;
							}
							if (nnnx > 0 && nl > 0 || nl == 0 && nnnx > 0) {
							}
							else if (nl > 0) {
								tc_apu2 = ty0_apu - tyu;
							}
							else if (first_y) {
								tc_apu2 = tc;
							}
							apu = (ty0_apu - tc_apu2) * L;

							elements[N2 + ii - nl - nnx] = (apu);
							temp += apu;
							ty0_apu -= tyu;
							apu_tempj -= ju;
						}
						apu_vary = false;
						nnn = 0;
						apu_varx2 = false;
						nnx = 0;
						first_y = false;
						first_x = false;
					}
					else {
						apu_tempj = tempj - ju;
						ty0_apu = ty0 - tyu;
						tc_apu = ty0_apu - tyu;
						for (int nl = static_cast<int>(nnn) - 1; nl >= 0; nl--) {
							indices[N2 + ii - nl] = apu_tempj*Nx + tempi + Nx * Ny * tempk;
							if (nl == 0)
								apu = (ty0_apu - tc) * L;
							else {
								apu = (ty0_apu - tc_apu) * L;
							}
							elements[N2 + ii - nl] = (apu);
							temp += apu;
							if (nnn > 1) {
								apu_tempj -= ju;
								ty0_apu -= tyu;
								tc_apu -= tyu;
							}
						}
						apu_vary = false;
						nnn = 0;
					}
				}
				else if (ii == Np - 1 && apu_varx && apu_vary) {
					bool apu_bool = false;
					if (final_x) {
						apu_tempi = tempi - iu;
						tx0_apu = tx0 - txu;
						apu_tempj = tempj;
						ty0_apu = ty0;
						if (tx0_apu - txu >= ty0_apu - tyu && nnn > 1) {
							tc_apu = tx0_apu - txu;
							apu_bool = true;
						}
						else if (nnn > 1)
							tc_apu = ty0_apu - tyu;
						else
							tc_apu = tc;
						apu = (tx0_apu - tc_apu) * L;
					}
					else {
						apu_tempj = tempj - ju;
						ty0_apu = ty0 - tyu;
						apu_tempi = tempi;
						tx0_apu = tx0;
						if (tx0_apu - txu >= ty0_apu - tyu && nnn > 1) {
							tc_apu = tx0_apu - txu;
							apu_bool = true;
						}
						else if (nnn > 1)
							tc_apu = ty0_apu - tyu;
						else
							tc_apu = tc;
						apu = (ty0_apu - tc_apu) * L;
					}
					for (int nl = static_cast<int>(nnn) - 1; nl >= 0; nl--) {
						
						indices[N2 + ii - nl] = apu_tempj*Nx + apu_tempi + Nx * Ny * tempk;
						elements[N2 + ii - nl] = (apu);
						temp += apu;
						if (apu_bool && nl > 0) {
							tx0_apu -= txu;
							apu_tempi -= iu;
							if (tx0_apu - txu >= ty0_apu - tyu && nl > 1) {
								tc_apu = tx0_apu - txu;
							}
							else if (nl > 1) {
								tc_apu = ty0_apu - tyu;
								apu_bool = false;
							}
							else
								tc_apu = tc;
							apu = (tx0_apu - tc_apu) * L;
						}
						else if (nl > 0) {
							ty0_apu -= tyu;
							apu_tempj -= ju;
							if (tx0_apu - txu >= ty0_apu - tyu && nl > 1) {
								tc_apu = tx0_apu - txu;
								apu_bool = true;
							}
							else if (nl > 1) {
								tc_apu = ty0_apu - tyu;
							}
							else
								tc_apu = tc;
							apu = (ty0_apu - tc_apu) * L;
						}
					}

				}
			}

			temp = 1. / temp;

			if (attenuation_correction) {

				double jelppi = 0.;

				for (unsigned int iii = 0; iii < Np; iii++) {
					jelppi += elements[N2 + iii] * -atten[indices[N2 + iii]];
				}
				temp = exp(jelppi) * temp;

			}

			for (unsigned int ii = 0; ii < Np; ii++) {
				elements[N2 + ii] = (fabs(elements[N2 + ii]) * temp);
			}
			lor[lo] = lor2[lo];
			continue;

		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray*prhs[])

{

	if (nrhs != 28)
		mexErrMsgTxt("Invalid number of input arguments.  There must be exactly 28.");

	if (nlhs != 1)
		mexErrMsgTxt("Invalid number of output arguments.  There must be exactly one.");


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

	vector<double> iij_vec(iij, iij + mxGetNumberOfElements(prhs[11]));

	vector<double> jjk_vec(jji, jji + mxGetNumberOfElements(prhs[12]));

	vector<double> kkj_vec(kkj, kkj + mxGetNumberOfElements(prhs[13]));

	vector<double> yy_vec(yy, yy + mxGetNumberOfElements(prhs[14]) - 1);

	vector<double> xx_vec(xx, xx + mxGetNumberOfElements(prhs[15]) - 1);

	int size_x = (int)mxGetScalar(prhs[16]);

	double zmax = (double)mxGetScalar(prhs[17]);

	double *atten = (double*)mxGetData(prhs[18]);

	unsigned int *lor2 = (unsigned int*)mxGetData(prhs[19]);

	int pituus = (int)mxGetScalar(prhs[20]);

	bool attenuation_correction = (bool)mxGetScalar(prhs[21]);

	uint16_t* lor1 = (uint16_t*)mxGetData(prhs[22]);

	uint64_t summa = (uint64_t)mxGetScalar(prhs[23]);

	uint16_t *L = (uint16_t*)mxGetData(prhs[24]);
	int numRows = (int)mxGetM(prhs[24]);

	int *pseudos = (int*)mxGetData(prhs[25]);
	int pRows = (int)mxGetM(prhs[25]);

	int det_per_ring = (int)mxGetScalar(prhs[26]);

	bool verbose = (bool)mxGetScalar(prhs[27]);

	size_t loop_var_par = pituus;

	double maxyy = *max_element(yy, yy + Ny + 1);
	double minyy = *min_element(yy, yy + Ny + 1);

	double maxxx = *max_element(xx, xx + Nx + 1);
	double minxx = *min_element(xx, xx + Nx + 1);

	mwSize N = Nx*Ny*Nz;
	mwSize nzmax = summa;
	mwSize rows = loop_var_par;

	plhs[0] = mxCreateSparse(N, rows, nzmax, mxREAL);
	double* elements = (double*)mxGetData(plhs[0]);
	mwIndex* indices = mxGetIr(plhs[0]);
	mwIndex* lor = mxGetJc(plhs[0]);

	clock_t time = clock();

	improved_siddon(loop_var_par, size_x, zmax, indices, elements, lor, maxyy, minyy, maxxx, minxx, xx_vec, iij_vec, jjk_vec, kkj_vec, yy_vec, atten, x, y, 
		z_det, Nx, Ny, Nz, d, dz, bx, by, bz, attenuation_correction, lor1, lor2, L, pseudos, pRows, det_per_ring);

	time = clock() - time;

	if (verbose) {
		mexPrintf("Improved Siddon took %f seconds\n", ((float)time) / CLOCKS_PER_SEC);
		mexEvalString("pause(.001);");
	}

}
