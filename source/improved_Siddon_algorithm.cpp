#include "mex.h"
#include <vector>
#include <algorithm>
#include <numeric>
#include <time.h>


using namespace std;

// this function was taken from: https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

	// initialize original index locations
	vector<size_t> idx(v.size());
	iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

	return idx;
}

int improved_siddon(const int loop_var_par, const int size_x, const double zmax, const int TotSinos, vector<int>& indices, vector<float>& elements, 
	vector<vector<int>>& lor, const double maxyy, const double minyy, const double maxxx, const double minxx, const vector<double>& xx_vec, 
	const vector<double>& z_det_vec, const vector<double>& iij_vec, const vector<double>& jjk_vec, const vector<double>& kkj_vec, const vector<double>& yy_vec, 
	const double* atten, const double* x, const double* y, const double* z_det, const int NSlices, const int Nx, const int Ny, const int Nz, const double d, 
	const double dz, const double bx, const double by, const double bz, const int pituus, const unsigned int *index, const bool attenuation_correction) {

	int ll = -1;
	int lz = -1;
	int lj = 0;
	int oo = -1;
	int lh = -1;

	const double bxf = bx + iij_vec.front() * d;
	const double byf = by + jjk_vec.front() * d;
	const double bxb = bx + iij_vec.back() * d;
	const double byb = by + jjk_vec.back() * d;
	const double bzf = bz + kkj_vec.front() * dz;
	const double bzb = bz + kkj_vec.back() * dz;

	for (int lo = 1; lo <= loop_var_par; lo++) {

		//mexPrintf("lo = %d\n", lo);

		//aa++;
		if ((lo - 1) % size_x == 0) {
			ll = -1;
			lz++;
		}

		//mexPrintf("aa = %d\n", aa);

		ll++;

		if (pituus > 1) {
			//mexPrintf("lo = %d\n", lo);
			if (lo != index[oo + 1]) {
				continue;
			}
		}
			

		oo++;

		//if (oo >= 7170302) {
		////if (oo % 100000 == 0) {
		//	mexPrintf("oo = %d\n", oo);
		//	mexPrintf("lo = %d\n", lo);
		//	mexEvalString("pause(.001);");
		//}
		//else if (oo > pituus - 1) {
			//mexPrintf("lo = %d\n", lo);

			//mexEvalString("pause(.001);");
		//}

		//}



// 		if (discard[lo - 1] == false)
// 			continue;

		double xs = x[ll];
		double xd = x[ll + size_x];
		double ys = y[ll];
		double yd = y[ll + size_x];
		double zs = z_det[lz];
		double zd = z_det[lz + TotSinos];
		

		double y_diff = (yd - ys);
		double x_diff = (xd - xs);
		double z_diff = (zd - zs);

		//mexPrintf("zs = %f\n", zs);
		//if (lo == 68) {
		//	mexPrintf("yd = %f\n", yd);
		//	mexPrintf("ys = %f\n", ys);
		//	mexPrintf("xd = %f\n", xd);
		//	mexPrintf("xs = %f\n", xs);
		//	//mexPrintf("zd = %f\n", z_diff);
		//}
		

		if (fabs(z_diff) < 1e-8) {

			int z_loop = static_cast<int>((zs / zmax)*(static_cast<double>(NSlices) - 1.));

			//if (lz == 0) {
			//	mexPrintf("z_loop = %d\n", z_loop);
			//	//mexPrintf("zd = %f\n", zd);
			//}

			//mexPrintf("erotus = %f\n", fabs(y_diff));
			//mexPrintf("lo2 = %d\n", lo);

			if (fabs(y_diff) < 1e-8) {

				//mexPrintf("lo = %d\n", lo);

				if (yd <= maxyy && yd >= minyy) {
					double minvalue = maxyy * 100.;
					int apu;
					for (int ii = 0; ii < Ny; ii++) {
						double temp = fabs(yy_vec[ii] - yd);
						if (temp < minvalue) {
							minvalue = temp;
							apu = ii;
						}
					}
					//vector<double> templ_ijk(pixelsize, d);
					double templ_ijk = d;
					vector<int> tempk(Nx, apu * Ny + z_loop * Nx * Ny);
					double temp = d * static_cast<double>(Nx);
					//if (lo == 1315)
					//	mexPrintf("pixelsize = %d\n", pixelsize);
					//for_each(templ_ijk.begin(), templ_ijk.end(), [&](int n) {
					//	temp += n;
					//});
					temp = 1. / temp;

					//if (lo == 1315)
					//	mexPrintf("temp = %f\n", temp);

					if (attenuation_correction) {

						double jelppi = 0.;

						for (int iii = 0; iii < tempk.size(); iii++) {
							jelppi += templ_ijk * -atten[tempk[iii] + iii];
							//if (lo == 7259089)
								//mexPrintf("atten = %f\n", -atten[tempk[iii] + iii]);
							//if (lo == 1315)
								//mexPrintf("tempk = %d\n", tempk[iii] + iii);
						}
						//if (lo == 7259089)
							//mexPrintf("jelppi = %f\n", jelppi);
						temp = exp(jelppi) * temp;

					}

					//if (lo == 1315)
					//	mexPrintf("jelppi = %f\n", jelppi);

					float element = static_cast<float>(templ_ijk * temp);

					for (int ii = 0; ii < tempk.size(); ii++) {
						indices.emplace_back((tempk[ii] + ii));
						lh++;
						elements.emplace_back(element);
						//if (elements.back() > 1.)
							//mexPrintf("lo = %d\n", lo);
					}
					lj++;

					lor[oo][0] = oo + 1;
					lor[oo][1] = Nx;
					//mexPrintf("lor = %d\n", lor[oo][0]);
					continue;
				}
				//mexPrintf("lo1 = %d\n", lo);
				continue;
			}
			else if (fabs(x_diff) < 1e-8) {

				//mexPrintf("lo = %d\n", lo);
				if (xd <= maxxx && xd >= minxx) {
					double minvalue = maxxx * 100.;
					int apu;
					for (int ii = 0; ii < Nx; ii++) {
						double temp = fabs(xx_vec[ii] - xd);
						if (temp < minvalue) {
							minvalue = temp;
							apu = ii;
						}
					}
					//vector<double> templ_ijk(pixelsize, d);
					double templ_ijk = d;
					vector<int> tempk(Ny, apu + z_loop * Nx * Ny);
					//double temp;
					double temp = d * static_cast<double>(Ny);
					//for_each(templ_ijk.begin(), templ_ijk.end(), [&](int n) {
					//	temp += n;
					//});
					temp = 1. / temp;

					if (attenuation_correction) {

						double jelppi = 0.;

						for (int iii = 0; iii < tempk.size(); iii++) {
							jelppi += templ_ijk * -atten[tempk[iii] + iii * Ny];
							//if (lo == 7259089)
								//mexPrintf("atten = %f\n", -atten[tempk[iii] + iii]);
							//if (lo == 7259089)
								//mexPrintf("templ_ijk = %f\n", templ_ijk[iii]);
						}
						temp = exp(jelppi) * temp;

					}
					//if (lo == 7259089)
						//mexPrintf("jelppi = %f\n", jelppi);

					float element = static_cast<float>(templ_ijk * temp);

					for (int ii = 0; ii < tempk.size(); ii++) {
						indices.emplace_back((tempk[ii] + ii * Ny));
						elements.emplace_back(element);
						lh++;
						//if (elements.back() > 1.)
							//mexPrintf("lo = %d\n", lo);
					}
					lj++;

					lor[oo][0] = oo + 1;
					lor[oo][1] = Ny;
					//mexPrintf("lor = %d\n", lor[oo][0]);
					continue;
				}
				//mexPrintf("lo2 = %d\n", lo);
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

			//mexPrintf("lo = %d\n", lo);

			if (tmin >= tmax) {
				//mexPrintf("tmin = %f\n", tmin);
				//mexPrintf("tmax = %f\n", tmax);
				//mexPrintf("lo3 = %d\n", lo);
				continue;
			}

			//if (lo == 342)
			//	mexPrintf("tmin = %f\n", tmin);

			//if (lo == 342)
			//	mexPrintf("tmax = %f\n", tmax);

			int imin, imax, jmin, jmax;
			double pxt, pyt;

			//vector<double> tx_n;
			//vector<double> ty_n;

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
				//apu = (imin:1 : imax) + 1;
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
				//vector<double> tx_n(tx.begin() + imin, tx.begin() + imax);
			}

			//if (lo == 342)
			//	mexPrintf("imin = %d\n", imin);

			//if (lo == 342)
			//	mexPrintf("imax = %d\n", imax);
			//mexPrintf("lo = %f\n", pxt);
			//mexPrintf("tmin = %f\n", tmin);
			//mexPrintf("tmax = %f\n", tmax);
			//if (imax > imin) {
			//	tx_n.reserve(imax + imin);
			//	tx_n.insert(tx_n.end(), tx.begin() + imin, tx.begin() + imax);
			//}
			//else if (imax == imin)
			//	tx_n.emplace_back(tx[imax]);

			//mexPrintf("lo2 = %d\n", lo);

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

			//if (lo == 342)
			//	mexPrintf("jmin = %d\n", jmin);

			//if (lo == 342)
			//	mexPrintf("jmax = %d\n", jmax);

			int Np = (imax - imin + 1) + (jmax - jmin + 1);

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

			vector<double> templ_ijk;
			vector<int> temp_koko;
			temp_koko.reserve(Np);
			templ_ijk.reserve(Np);

			for (int ii = 0; ii < Np; ii++) {

				temp_koko.emplace_back(tempj*Nx + tempi + tempk);

				if (tx0 < ty0) {

					
					templ_ijk.emplace_back((tx0 - tc) * L);

					//if (ii == 0) {
					//	indices[lo][ii] = (tempj + tempi + tempk);
					//	templ_ijk[ii] = (tx0 - tc) * L;
					//}
					//else {
					//	indices[lo].emplace_back((tempj + tempi + tempk));
					//	templ_ijk.emplace_back((tx0 - tc) * L);
					//}

					tempi += iu;
					tc = tx0;
					tx0 += txu;

					temp += templ_ijk[ii];

					//mexPrintf("ii = %d\n", (tempj * Ny + tempi + tempk * Ny * Nx));


				}
				else {

					templ_ijk.emplace_back((ty0 - tc) * L);

					tempj += ju;
					//tempj *= Nx;
					tc = ty0;
					ty0 += tyu;

					temp += templ_ijk[ii];
				}
				//mexPrintf("temp-koko = %d\n", temp_koko[ii]);
				//if (lo == 7319201) {
				//	mexPrintf("tempk = %d\n", tempk);
				//	mexPrintf("tempj = %d\n", tempj);
				//	mexPrintf("tempi = %d\n", tempi);
				//}
				if (tempj < 0 || tempi < 0 || tempi >= Nx || tempj >= Ny)
					break;


			}

			//mexPrintf("lo4 = %d\n", lo);

			temp = 1. / temp;

			//if (lo == 68)
			//	mexPrintf("temp = %f\n", temp);

			if (attenuation_correction) {

				double jelppi = 0.;

				for (int iii = 0; iii < templ_ijk.size(); iii++) {
					//if (lo == 7319201)
					//	mexPrintf("temp_koko = %d\n", temp_koko[iii]);
					jelppi += templ_ijk[iii] * -atten[temp_koko[iii]];
					//jelppi += temp_all.at(iii).templ_ijk.front() * -atten[temp_all.at(iii).temp_koko.front()];
				}
				temp = exp(jelppi) * temp;

			}

			auto i = sort_indexes(temp_koko);

			//mexPrintf("lo5 = %d\n", lo);

			for (int ii = 0; ii < templ_ijk.size(); ii++) {
				//if (lo == 121)
				//	mexPrintf("element = %f\n", templ_ijk[i[ii]]);
				elements.emplace_back(static_cast<float>(templ_ijk[i[ii]] * temp));
				indices.emplace_back(temp_koko[i[ii]]);
				lh++;
				//if (temp_koko[i[ii]] < (-(int)2147483647))
				//	mexPrintf("lo5 = %d\n", lo);
				//if (lo == 342) {
				//	mexPrintf("temp_koko = %d\n", temp_koko[i[ii]]);
				//	mexPrintf("templ_ijk = %f\n", templ_ijk[i[ii]]);
				//}
				
				//if (elements.back() > 1.)
				//	mexPrintf("lo = %d\n", lo);
			}

			//if (indices[lh - 3] == 16127 && indices[lh - 2] == 16255 && indices[lh - 1] == 16382 && indices[lh] == 16383)
			//	mexPrintf("lo5 = %d\n", lo);


			lor[oo][0] = oo + 1;
			lor[oo][1] = static_cast<int>(templ_ijk.size());
			lj++;
			//mexPrintf("lo = %d\n", lor[oo][0]);
			//mexPrintf("lo6 = %d\n", lo);
			continue;
		}
		else {

			//if (lo > 41730900)
				//mexPrintf("lo1 = %d\n", lo);
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
						//mexPrintf("tmin = %f\n", tmin);
						//mexPrintf("tmax = %f\n", tmax);
						//mexPrintf("lo4 = %d\n", lo);
						continue;
					}

					int imin, imax, kmin, kmax;
					double pxt, pzt;

					int iu, ku;

					//vector<double> tx_n;
					//vector<double> tz_n;

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

					double L = sqrt((x_diff*x_diff + z_diff*z_diff));

					int tempi, tempj, tempk;
					double apu2, apu1;

					double pt = ((min(tx0, tz0) + tmin) / 2.);

					tempi = static_cast<int>(floor(((xs + pt * x_diff) - bx) / d));
					tempk = static_cast<int>(floor(((zs + pt * z_diff) - bz) / dz));

					double txu = d / fabs(x_diff);
					double tzu = dz / fabs(z_diff);

					vector<double> templ_ijk;
					vector<int> temp_koko;

					temp_koko.reserve(Np);
					templ_ijk.reserve(Np);

					for (int ii = 0; ii < Ny; ii++) {
						apu1 = fabs(yy_vec[ii] - yd);
						if (ii > 0 && apu1 < apu2 || ii == 0) {
							tempj = ii;
							apu2 = apu1;
						}

					}

					tempj = tempj * Nx;

					double tc = tmin;

					double temp = 0.;

					for (int ii = 0; ii < Np; ii++) {

						//mexPrintf("ii = %d\n", (tempj * Ny + tempi + tempk * Ny * Nx));

						if (tx0 < tz0) {

							temp_koko.emplace_back(tempj + tempi + Nx * Ny * tempk);
							templ_ijk.emplace_back((tx0 - tc) * L);


							tempi += iu;
							tc = tx0;
							tx0 += txu;

							temp += templ_ijk[ii];
						}
						else {

							temp_koko.emplace_back(tempj + tempi + Nx * Ny * tempk);
							templ_ijk.emplace_back((tz0 - tc) * L);

							tempk += ku;
							//tempk *= (Nx * Ny);
							tc = tz0;
							tz0 += tzu;

							temp += templ_ijk[ii];
						}

						if (tempk < 0 || tempi < 0 || tempi >= Nx || tempk >= Nz)
							break;

					}

					temp = 1. / temp;

					if (attenuation_correction) {

						double jelppi = 0.;

						for (int iii = 0; iii < templ_ijk.size(); iii++) {
							//if (lo == 7319201)
							//	mexPrintf("temp_koko = %d\n", temp_koko[iii]);
							jelppi += templ_ijk[iii] * -atten[temp_koko[iii]];
							//jelppi += temp_all.at(iii).templ_ijk.front() * -atten[temp_all.at(iii).temp_koko.front()];
						}
						temp = exp(jelppi) * temp;

					}

					auto i = sort_indexes(temp_koko);

					//mexPrintf("lo5 = %d\n", lo);

					for (int ii = 0; ii < templ_ijk.size(); ii++) {
						elements.emplace_back(static_cast<float>(templ_ijk[i[ii]] * temp));
						indices.emplace_back(temp_koko[i[ii]]);
						lh++;
						//if (temp_koko[i[ii]] < (-(int)2147483647))
							//mexPrintf("lo5 = %d\n", lo);
					}


					lor[oo][0] = oo + 1;
					lor[oo][1] = static_cast<int>(templ_ijk.size());
					lj++;

					continue;

				}
				//mexPrintf("lo5 = %d\n", lo);
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
						//mexPrintf("tmin = %f\n", tmin);
						//mexPrintf("tmax = %f\n", tmax);
						//mexPrintf("lo6 = %d\n", lo);
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

					vector<double> templ_ijk;
					vector<int> temp_koko;

					temp_koko.reserve(Np);
					templ_ijk.reserve(Np);

					for (int ii = 0; ii < Nx; ii++) {
						apu1 = fabs(xx_vec[ii] - xd);
						if (ii > 0 && apu1 < apu2 || ii == 0) {
							tempi = ii;
							apu2 = apu1;
						}

					}

					//tempi = tempi * Ny;

					for (int ii = 0; ii < Np; ii++) {


						if (tz0 < ty0) {

							temp_koko.emplace_back(tempj*Nx + tempi + Nx * Ny * tempk);
							templ_ijk.emplace_back((tz0 - tc) * L);

							tempk += ku;
							//tempk *= (Nx * Ny);
							tc = tz0;
							tz0 += tzu;

							temp = temp + templ_ijk[ii];

						}
						else {

							temp_koko.emplace_back(tempj*Nx + tempi + Nx * Ny * tempk);
							templ_ijk.emplace_back((ty0 - tc) * L);

							tempj += ju;
							//tempj *= Nx;
							tc = ty0;
							ty0 += tyu;

							temp += templ_ijk[ii];
						}
						if (tempj < 0 || tempk < 0 || tempk >= Nz || tempj >= Ny)
							break;

					}

					temp = 1. / temp;

					if (attenuation_correction) {

						double jelppi = 0.;

						for (int iii = 0; iii < templ_ijk.size(); iii++) {
							//if (lo == 7319201)
							//	mexPrintf("temp_koko = %d\n", temp_koko[iii]);
							jelppi += templ_ijk[iii] * -atten[temp_koko[iii]];
							//jelppi += temp_all.at(iii).templ_ijk.front() * -atten[temp_all.at(iii).temp_koko.front()];
						}
						temp = exp(jelppi) * temp;

					}

					auto i = sort_indexes(temp_koko);

					//mexPrintf("lo5 = %d\n", lo);

					for (int ii = 0; ii < templ_ijk.size(); ii++) {
						elements.emplace_back(static_cast<float>(templ_ijk[i[ii]] * temp));
						indices.emplace_back(temp_koko[i[ii]]);
						lh++;
						//if (temp_koko[i[ii]] < (-(int)2147483647))
							//mexPrintf("lo5 = %d\n", lo);
					}


					lor[oo][0] = oo + 1;
					lor[oo][1] = static_cast<int>(templ_ijk.size());
					lj++;

					continue;

				}
				//mexPrintf("lo7 = %d\n", lo);
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
				//mexPrintf("tmin = %f\n", tmin);
				//mexPrintf("tmax = %f\n", tmax);
				//mexPrintf("lo8 = %d\n", lo);
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

			/*if (lo == 87693440) {
				mexPrintf("z_diff = %f\n", z_diff);
				mexPrintf("kmax = %d\n", kmax);
				mexPrintf("dz = %f\n", dz);
				mexPrintf("zs = %f\n", zs);
				mexPrintf("zd = %f\n", zd);
				mexPrintf("bz = %f\n", bz);
			}*/

			int Np = (kmax - kmin + 1) + (jmax - jmin + 1) + (imax - imin + 1);

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

			vector<double> templ_ijk;
			vector<int> temp_koko;

			temp_koko.reserve(Np);
			templ_ijk.reserve(Np);

			//if (lo == 29668657) {
			//	mexPrintf("tz0 = %f\n", tz0);
			//	mexPrintf("z_diff = %f\n", z_diff);
			//	mexPrintf("zs = %f\n", zs);
			//	mexPrintf("bz = %f\n", bz);
			//	mexPrintf("dz = %f\n", dz);
			//	mexPrintf("kmax = %d\n", kmax);
			//	mexPrintf("kmin = %d\n", kmin);
			//}


			for (int ii = 0; ii < Np; ii++) {


				if (tz0 < ty0 && tz0 < tx0) {

					temp_koko.emplace_back(tempj*Nx + tempi + Nx * Ny * tempk);
					templ_ijk.emplace_back((tz0 - tc) * L);

					tempk += ku;
					//tempk *= (Nx * Ny);
					tc = tz0;
					tz0 += tzu;

					temp += templ_ijk[ii];
				}
				else if (ty0 < tx0) {

					temp_koko.emplace_back(tempj*Nx + tempi + Nx * Ny * tempk);
					templ_ijk.emplace_back((ty0 - tc) * L);

					tempj += ju;
					//tempj *= Nx;
					tc = ty0;
					ty0 += tyu;

					temp += templ_ijk[ii];
				}
				else {

					temp_koko.emplace_back(tempj*Nx + tempi + Nx * Ny * tempk);
					templ_ijk.emplace_back((tx0 - tc) * L);

					tempi += iu;
					tc = tx0;
					tx0 += txu;

					temp += templ_ijk[ii];
				}
				//if (lo == 87693440) {
				//	mexPrintf("tempj = %d\n", tempj);
				//	mexPrintf("tempi= %d\n", tempi);
				//	mexPrintf("tempk = %d\n", tempk);
				//}
				if (tempj < 0 || tempi < 0 || tempk < 0 || tempi >= Nx || tempj >= Ny || tempk >= Nz)
					break;

			}

			

			temp = 1. / temp;

			if (attenuation_correction) {

				double jelppi = 0.;

				for (int iii = 0; iii < templ_ijk.size(); iii++) {
					//if (lo == 7319201)
					//	mexPrintf("temp_koko = %d\n", temp_koko[iii]);
					jelppi += templ_ijk[iii] * -atten[temp_koko[iii]];
					//jelppi += temp_all.at(iii).templ_ijk.front() * -atten[temp_all.at(iii).temp_koko.front()];
				}
				temp = exp(jelppi) * temp;

			}

			auto i = sort_indexes(temp_koko);

			//mexPrintf("lo5 = %d\n", lo);

			for (int ii = 0; ii < templ_ijk.size(); ii++) {
				elements.emplace_back(static_cast<float>(templ_ijk[i[ii]] * temp));
				indices.emplace_back(temp_koko[i[ii]]);
				lh++;
				//if (temp_koko[i[ii]] < (-(int)2147483647))
					//mexPrintf("lo5 = %d\n", lo);
			}

			//if (indices[lh - 3] == 382211 && indices[lh -2] == 386176 && indices[lh - 1] == 386240 && indices[lh] == 386241)
			//	mexPrintf("lo5 = %d\n", lo);


			lor[oo][0] = oo + 1;
			lor[oo][1] = static_cast<int>(templ_ijk.size());
			lj++;

			continue;
		}
		//mexPrintf("lo9 = %d\n", lo);
	}
	return lj;
}

void mexFunction(int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray*prhs[])

{
	if (nrhs != 27)
		mexErrMsgTxt("Invalid number of input arguments.  There must be exactly 27.");

	if (nlhs != 3)
		mexErrMsgTxt("Invalid number of output arguments.  There must be exactly three.");

	
	bool verbose = (bool)mxGetScalar(prhs[0]);

	int Ny = (int)mxGetScalar(prhs[1]);

	int Nx = (int)mxGetScalar(prhs[2]);

	int Nz = (int)mxGetScalar(prhs[3]);

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

	int NSinos = (int)mxGetScalar(prhs[17]);

	int NSlices = (int)mxGetScalar(prhs[18]);

	//int NSlices = (int)mxGetScalar(prhs[18]);

	vector<double> z_det_vec(z_det, z_det + mxGetNumberOfElements(prhs[9]));

	//vector<double> x_vec(x, x + mxGetNumberOfElements(prhs[10]));

	//vector<double> y_vec(y, y + mxGetNumberOfElements(prhs[11]));

	vector<double> iij_vec(iij, iij + mxGetNumberOfElements(prhs[12]));

	//for (int kk = 0; kk < mxGetNumberOfElements(prhs[13]); kk++)
		//mexPrintf("jji = %f\n", jji[kk]);

	vector<double> jjk_vec(jji, jji + mxGetNumberOfElements(prhs[13]));

	vector<double> kkj_vec(kkj, kkj + mxGetNumberOfElements(prhs[14]));

	vector<double> yy_vec(yy, yy + mxGetNumberOfElements(prhs[15]) - 1);

	vector<double> xx_vec(xx, xx + mxGetNumberOfElements(prhs[16]) - 1);

	int size_x = (int)mxGetScalar(prhs[19]);

	double zmax = (double)mxGetScalar(prhs[20]);

	int TotSinos = (int)mxGetScalar(prhs[21]);

	int ind_size = (int)mxGetScalar(prhs[22]);

	double *atten = (double*)mxGetData(prhs[23]);

	unsigned int *index = (unsigned int*)mxGetData(prhs[24]);
	int index_size = mxGetNumberOfElements(prhs[24]);

	int pituus = (int)mxGetScalar(prhs[25]);

	bool attenuation_correction = (bool)mxGetScalar(prhs[26]);

// 	bool *discard = (bool*)mxGetData(prhs[24]);

	int loop_var_par = 1;

	//mexPrintf("pituus = %d\n", pituus);

	//mexEvalString("pause(.001);");

	vector<vector<int>> lor;

	if (index_size > 1) {
		loop_var_par = index[pituus - 1];
		lor.assign(pituus, vector<int>(2, 0));
	}
	else {
		loop_var_par = NSinos * size_x;
		lor.assign(loop_var_par, vector<int>(2, 0));
	}

	//mexPrintf("loop_var_par = %d\n", loop_var_par);

	//mexEvalString("pause(.001);");

	//mexPrintf("lor.size() = %d\n", lor.size());

	//mexEvalString("pause(.001);");

	//mexPrintf("lor[0].size() = %d\n", lor[0].size());

	//mexEvalString("pause(.001);");

	vector<int> indices;
	//vector<vector<int>> indices;

	//indices.resize(loop_var_par);

	vector<float> elements;

	indices.reserve(ind_size);
	elements.reserve(ind_size);


	
	//vector<vector<int>> lor;

	double maxyy = *max_element(yy, yy + Ny + 1);
	double minyy = *min_element(yy, yy + Ny + 1);

	double maxxx = *max_element(xx, xx + Nx + 1);
	double minxx = *min_element(xx, xx + Nx + 1);

	//mexPrintf("max = %f\n", maxyy);

	//mexEvalString("pause(.001);");
	//mexPrintf("min = %f\n", minyy);

	clock_t time = clock();

	int lj = improved_siddon(loop_var_par, size_x, zmax, TotSinos, indices, elements, lor, maxyy, minyy, maxxx, minxx, xx_vec, z_det_vec, iij_vec, jjk_vec, 
		kkj_vec, yy_vec, atten, x, y, z_det, NSlices, Nx, Ny, Nz, d, dz, bx, by, bz, index_size, index, attenuation_correction);
    
	time = clock() - time;

	if (verbose) {
		mexPrintf("Function elapsed time is %f seconds\n", ((float)time) / CLOCKS_PER_SEC);
		mexEvalString("pause(.001);");
	}
	//vector<int> global_vector;

	//mexPrintf("lo = %d\n", lor[2][0]);
	//mexPrintf("lo = %d\n", lor[loop_var_par - 1][1]);

	//mexPrintf("lj = %d\n", lj);

	//mexEvalString("pause(.001);");

	size_t outSize1 = lj * 2;
	size_t outSize2 = 1;

	//plhs[0] = mxCreateNumericMatrix(outSize1, outSize2, mxINT32_CLASS, mxREAL);
	plhs[0] = mxCreateNumericMatrix(lor.size() * 2, outSize2, mxINT32_CLASS, mxREAL);

	int* outputMatrix = (int *)mxGetData(plhs[0]);

	plhs[1] = mxCreateNumericMatrix(indices.size(), outSize2, mxINT32_CLASS, mxREAL);

	int* outputMatrix2 = (int *)mxGetData(plhs[1]);

	plhs[2] = mxCreateNumericMatrix(elements.size(), outSize2, mxSINGLE_CLASS, mxREAL);

	float* outputMatrix3 = (float *)mxGetData(plhs[2]);

	//plhs[3] = mxCreateLogicalMatrix(loop_var_par, outSize2);

	//bool* discard = (bool*)mxGetLogicals(plhs[3]);

	//plhs[0] = mxCreateNumericMatrix(outSize1, outSize2, mxSINGLE_CLASS, mxREAL);

	//double* outputMatrix = (double *)mxGetData(plhs[0]);

	//int* outputMatrix = lor.data();

	//mexPrintf("lor.size() = %d\n", lor.size());

	//mexEvalString("pause(.001);");

	time = clock();


	int lb = 0;
	for (int col = 0; col < 2; col++) {
		int la = col*lor.size();
		//int il = col*lor.size();
		for (int row = 0; row < lor.size(); row++) {
			//if (lor[row][0] == 0) {
			//	lb++;
			//	continue;
			//}
			//if (col == 0)
			//	lor[row][0] -= lb;
			//outputMatrix[la + col*lj] = lor[row][col];
			outputMatrix[la] = lor[row][col];
			la++;
		}
	}

	time = clock() - time;

	if (verbose)
		mexPrintf("lor copy elapsed time is %f seconds\n", ((float)time) / CLOCKS_PER_SEC);

	time = clock();

	//lor.~vector();
	lor.erase(lor.begin(), lor.end());
	lor.shrink_to_fit();

	time = clock() - time;
	if (verbose)
		mexPrintf("lor erase elapsed time is %f seconds\n", ((float)time) / CLOCKS_PER_SEC);

	//outputMatrix2 = indices.data();
	//outputMatrix3 = elements.data();

	//std::memcpy(mxGetData(plhs[1]), &indices[0], (indices.size())*(outSize2) * sizeof(int));

	time = clock();

	copy(indices.begin(), indices.end(), outputMatrix2);
	indices.erase(indices.begin(), indices.end());
	indices.shrink_to_fit();

	time = clock() - time;
	if (verbose)
		mexPrintf("indices copy and erase elapsed time is %f seconds\n", ((float)time) / CLOCKS_PER_SEC);
	//indices.~vector();
	time = clock();

	copy(elements.begin(), elements.end(), outputMatrix3);

	time = clock() - time;
	if (verbose)
		mexPrintf("elements copy elapsed time is %f seconds\n", ((float)time) / CLOCKS_PER_SEC);

	//for (int row = 0; row < indices.size(); row++) {
	//	outputMatrix2[row] = indices[row];
	//	outputMatrix3[row] = elements[row];
	//}

	//for (auto & buffer : lor) {
	//	move(buffer.begin(), buffer.end(), back_inserter(global_vector));
	//}

	

	//std::memcpy(mxGetData(plhs[0]), &lor[0][0], (outSize1)*(outSize2) * sizeof(int));
}