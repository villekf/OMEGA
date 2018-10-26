#include "mex.h"
#include <vector>
#include <algorithm>


using namespace std;

void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray*prhs[])

{
	if (nrhs != 21)
		mexErrMsgTxt("Invalid number of input arguments.  There must be exactly 21.");

	if (nlhs != 2)
		mexErrMsgTxt("Invalid number of output arguments.  There must be exactly two.");

	
	int Nx = (int)mxGetScalar(prhs[0]);

	int block1 = (int)mxGetScalar(prhs[1]);

	int blocks = (int)mxGetScalar(prhs[2]);

	int det_per_ring = (int)mxGetScalar(prhs[3]);

	double d = (double)mxGetScalar(prhs[4]);

	double dz = (double)mxGetScalar(prhs[5]);

	double by = (double)mxGetScalar(prhs[6]);

	double bx = (double)mxGetScalar(prhs[7]);

	double bz = (double)mxGetScalar(prhs[8]);

	double *z_det = (double*)mxGetData(prhs[9]);

	double *x = (double*)mxGetData(prhs[10]);

	double *y = (double*)mxGetData(prhs[11]);

	double *iij = (double*)mxGetData(prhs[12]);

	double *jji = (double*)mxGetData(prhs[13]);

	double *kkj = (double*)mxGetData(prhs[14]);

	double *yy = (double*)mxGetData(prhs[15]);

	double *xx = (double*)mxGetData(prhs[16]);

	uint16_t *L = (uint16_t*)mxGetData(prhs[17]);
	int numRows = (int)mxGetM(prhs[17]);

	int Ny = (int)mxGetScalar(prhs[18]);

	int Nz = (int)mxGetScalar(prhs[19]);

	int *pseudos = (int*)mxGetData(prhs[20]);
	int pRows = (int)mxGetM(prhs[20]);

	const vector<double> z_det_vec(z_det, z_det + mxGetNumberOfElements(prhs[9]));

	const vector<double> iij_vec(iij, iij + mxGetNumberOfElements(prhs[12]));


	const vector<double> jjk_vec(jji, jji + mxGetNumberOfElements(prhs[13]));

	const vector<double> kkj_vec(kkj, kkj + mxGetNumberOfElements(prhs[14]));

	const vector<double> yy_vec(yy, yy + mxGetNumberOfElements(prhs[15]) - 1);

	const vector<double> xx_vec(xx, xx + mxGetNumberOfElements(prhs[16]) - 1);

	const int loop_var_par = numRows/2;

	plhs[0] = mxCreateLogicalMatrix(loop_var_par, 1);

	bool* discard = (bool*)mxGetLogicals(plhs[0]);

	plhs[1] = mxCreateNumericMatrix(loop_var_par, 1, mxUINT16_CLASS, mxREAL);

	uint16_t* lor = (uint16_t*)mxGetData(plhs[1]);

	int ll = 0;

	const double maxyy = *max_element(yy, yy + Ny + 1);
	const double minyy = *min_element(yy, yy + Ny + 1);

	const double maxxx = *max_element(xx, xx + Nx + 1);
	const double minxx = *min_element(xx, xx + Nx + 1);

	int ps;
	double zs, zd;

	const double bxf = bx + iij_vec.front() * d;
	const double byf = by + jjk_vec.front() * d;
	const double bxb = bx + iij_vec.back() * d;
	const double byb = by + jjk_vec.back() * d;
	const double bzf = bz + kkj_vec.front() * dz;
	const double bzb = bz + kkj_vec.back() * dz;
	
	for (int lo = 1; lo <= loop_var_par; lo++) {


		int detektorit1 = static_cast<int>(L[ll++]);
		int detektorit2 = static_cast<int>(L[ll++]);


		if (detektorit1 == detektorit2)
			continue;

		int loop1 = 1 + ((detektorit1 - 1) / det_per_ring);
		int loop2 = 1 + ((detektorit2 - 1) / det_per_ring);


		

		if (loop1 == loop2) {
			if (loop1 > blocks || loop1 < block1 || loop2 > blocks || loop2 < block1)
				continue;
			zs = z_det_vec[loop1 - 1];
			zd = zs;
		}
		else {
			if (loop1 > blocks && loop2 > blocks || loop1 < block1 && loop2 < block1)
				continue;
			zs = z_det_vec[loop1 - 1];
			zd = z_det_vec[loop2 - 1];
		}

		

		ps = 0;

		for (int kk = 0; kk < pRows; kk++) {
			if (kk + 1 < pRows) {
				if (loop1 >= pseudos[kk] && loop1 < pseudos[kk + 1]) {
					zs = z_det_vec[loop1 + ps];
					break;
				}
				else
					ps++;
			}
			else {
				if (loop1 >= pseudos[kk])
					zs = z_det_vec[loop1 + ps];
			}
		}

		ps = 0;

		for (int kk = 0; kk < pRows; kk++) {
			if (kk + 1 < pRows) {
				if (loop2 >= pseudos[kk] && loop2 < pseudos[kk + 1]) {
					zd = z_det_vec[loop2 + ps];
					break;
				}
				else
					ps++;
			}
			else {
				if (loop2 >= pseudos[kk])
					zd = z_det_vec[loop2 + ps];
			}
		}

		double xs = x[detektorit1 - det_per_ring*(loop1 - 1) - 1];
		double xd = x[detektorit2 - det_per_ring*(loop2 - 1) - 1];
		double ys = y[detektorit1 - det_per_ring*(loop1 - 1) - 1];
		double yd = y[detektorit2 - det_per_ring*(loop2 - 1) - 1];
		double x_diff = xd - xs;
		double y_diff = yd - ys;
		double z_diff = zd - zs;

		if (fabs(z_diff) < 1e-8) {


			if (fabs(y_diff) < 1e-8) {


				if (yd <= maxyy && yd >= minyy) {
					discard[lo - 1] = true;
					lor[lo - 1] = static_cast<uint16_t>(Nx);
					continue;
				}
				continue;
			}
			else if (fabs(x_diff) < 1e-8) {

				if (xd <= maxxx && xd >= minxx) {
					discard[lo - 1] = true;
					lor[lo - 1] = static_cast<uint16_t>(Ny);
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

			if (tmin >= tmax) {
				continue;
			}

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

			int Np = (imax - imin + 1) + (jmax - jmin + 1);

			int tempi, tempj;

			double pt = ((min(tx0, ty0) + tmin) / 2.);

			tempi = static_cast<int>(floor(((xs + pt * x_diff) - bx) / d));
			tempj = static_cast<int>(floor(((ys + pt * y_diff) - by) / d));

			double txu = d / fabs(x_diff);
			double tyu = d / fabs(y_diff);

			double tc = tmin;

			uint16_t temp_koko = 0;

			for (int ii = 0; ii < Np; ii++) {

				temp_koko++;

				if (tx0 < ty0) {

					tempi += iu;
					tc = tx0;
					tx0 += txu;

				}
				else {

					tempj += ju;
					tc = ty0;
					ty0 += tyu;
				}
				if (tempj < 0 || tempi < 0 || tempi >= Nx || tempj >= Ny)
					break;


			}

			discard[lo - 1] = true;
			lor[lo - 1] = temp_koko;

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

					if (tmin >= tmax) {
						continue;
					}

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

					int Np = (imax - imin + 1) + (kmax - kmin + 1);

					int tempi, tempj, tempk;
					double apu2, apu1;

					double pt = ((min(tx0, tz0) + tmin) / 2.);

					tempi = static_cast<int>(floor(((xs + pt * x_diff) - bx) / d));
					tempk = static_cast<int>(floor(((zs + pt * z_diff) - bz) / dz));

					double txu = d / fabs(x_diff);
					double tzu = dz / fabs(z_diff);

					uint16_t temp_koko = 0;


					for (int ii = 0; ii < Ny; ii++) {
						apu1 = fabs(yy_vec[ii] - yd);
						if (ii > 0 && apu1 < apu2 || ii == 0) {
							tempj = ii;
							apu2 = apu1;
						}

					}

					tempj = tempj * Nx;

					double tc = tmin;


					for (int ii = 0; ii < Np; ii++) {



						if (tx0 < tz0) {

							temp_koko++;

							tempi += iu;
							tc = tx0;
							tx0 += txu;

						}
						else {
							temp_koko++;


							tempk += ku;
							tc = tz0;
							tz0 += tzu;

						}

						if (tempk < 0 || tempi < 0 || tempi >= Nx || tempk >= Nz)
							break;

					}



					discard[lo - 1] = true;
					lor[lo - 1] = temp_koko;
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

					if (tmin >= tmax) {
						continue;
					}

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

					int Np = (kmax - kmin + 1) + (jmax - jmin + 1);


					int tempi, tempj, tempk;
					double apu2, apu1;

					double pt = ((min(tz0, ty0) + tmin) / 2.);

					tempk = static_cast<int>(floor(((zs + pt * z_diff) - bz) / dz));
					tempj = static_cast<int>(floor(((ys + pt * y_diff) - by) / d));

					double tzu = dz / fabs(z_diff);
					double tyu = d / fabs(y_diff);

					double tc = tmin;


					uint16_t temp_koko = 0;


					for (int ii = 0; ii < Nx; ii++) {
						apu1 = fabs(xx_vec[ii] - xd);
						if (ii > 0 && apu1 < apu2 || ii == 0) {
							tempi = ii;
							apu2 = apu1;
						}

					}


					for (int ii = 0; ii < Np; ii++) {


						if (tz0 < ty0) {

							temp_koko++;

							tempk += ku;
							tc = tz0;
							tz0 += tzu;

						}
						else {

							temp_koko++;

							tempj += ju;
							tc = ty0;
							ty0 += tyu;

						}
						if (tempj < 0 || tempk < 0 || tempk >= Nz || tempj >= Ny)
							break;

					}

					discard[lo - 1] = true;
					lor[lo - 1] = temp_koko;

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

			if (tmin >= tmax) {
				continue;
			}

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

			int Np = (kmax - kmin + 1) + (jmax - jmin + 1) + (imax - imin + 1);

			int tempi, tempj, tempk;

			double pt = ((min(min(tz0, ty0), tx0) + tmin) / 2.);

			tempk = static_cast<int>(floor(((zs + pt * z_diff) - bz) / dz));
			tempj = static_cast<int>(floor(((ys + pt * y_diff) - by) / d));
			tempi = static_cast<int>(floor(((xs + pt * x_diff) - bx) / d));


			double tzu = dz / fabs(z_diff);
			double tyu = d / fabs(y_diff);
			double txu = d / fabs(x_diff);

			double tc = tmin;
			uint16_t temp_koko = 0;

			for (int ii = 0; ii < Np; ii++) {


				if (tz0 < ty0 && tz0 < tx0) {

					temp_koko++;

					tempk += ku;
					tc = tz0;
					tz0 += tzu;
				}
				else if (ty0 < tx0) {

					temp_koko++;

					tempj += ju;
					tc = ty0;
					ty0 += tyu;
				}
				else {

					temp_koko++;

					tempi += iu;
					tc = tx0;
					tx0 += txu;
				}
				if (tempj < 0 || tempi < 0 || tempk < 0 || tempi >= Nx || tempj >= Ny || tempk >= Nz)
					break;

			}

			discard[lo - 1] = true;
			lor[lo - 1] = temp_koko;

			continue;

		}

	}
	
	return;
}