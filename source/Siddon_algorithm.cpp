#include "mex.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>


using namespace std;

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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
	if (nrhs != 27)
		mexErrMsgTxt("Invalid number of input arguments.  There must be exactly 27.");

	if (nlhs != 3)
		mexErrMsgTxt("Invalid number of output arguments.  There must be exactly three.");

	
	int pixelsize = (int)mxGetScalar(prhs[0]);

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

	int *index = (int*)mxGetData(prhs[24]);

	int pituus = (int)mxGetScalar(prhs[25]);

	bool attenuation_correction = (bool)mxGetScalar(prhs[26]);

	int loop_var_par;

	if (pituus > 1)
		loop_var_par = index[pituus - 1];
	else
		loop_var_par = NSinos * size_x;

	vector<int> indices;
	//vector<vector<int>> indices;

	//indices.resize(loop_var_par);

	vector<float> elements;

	indices.reserve(ind_size);
	elements.reserve(ind_size);

	vector<vector<int>> lor(loop_var_par, vector<int>(2, 0));
	//vector<vector<int>> lor;

	int ll = -1;
	int lz = -1;
	int lj = 0;
	int la = 0;
	int oo = 0;

	double maxyy = *max_element(yy, yy + Ny + 1);
	double minyy = *min_element(yy, yy + Ny + 1);

	double maxxx = *max_element(xx, xx + Nx + 1);
	double minxx = *min_element(xx, xx + Nx + 1);

	//mexPrintf("max = %f\n", maxyy);
	//mexPrintf("min = %f\n", minyy);

	

	for (int lo = 1; lo <= loop_var_par; lo++) {

		//if (lo > 3629023)
		//	mexPrintf("lo1 = %d\n", lo);

		if ((lo - 1) % size_x == 0) {
			ll = -1;
			lz++;
		}

		//if (lo > 3629023)
		//	mexPrintf("lo2 = %d\n", lo);

		//mexPrintf("z_loop = %d\n", lz);

		ll++;

		if (pituus > 1) {
			if (lo = !(index[oo] - 1))
				continue;
			else
				oo++;
		}

		double xs = x[ll];
		double xd = x[ll + size_x];
		double ys = y[ll];
		double yd = y[ll + size_x];
		double zs = z_det[lz];
		double zd = z_det[lz + TotSinos];
		int z_loop = static_cast<int>((zs / zmax)*(NSlices - 1) + 0.5);

		//if (lo > 3629023)
		//	mexPrintf("lo3 = %d\n", lo);

		//if (lo == 1008289) {
		//	mexPrintf("z_loop = %d\n", z_loop);
		//	mexPrintf("zs = %d\n", (int)((zs / zmax)*(NSlices - 1)+ 0.5));
		//	mexPrintf("zmax = %f\n", zmax);
		//}

		double y_diff = (yd - ys);
		double x_diff = (xd - xs);
		double z_diff = (zd - zs);

		//if (lo > 3629023)
		//	mexPrintf("lo4 = %d\n", lo);

		//mexPrintf("zs = %f\n", zs);
		//mexPrintf("zd = %f\n", z_diff);
		//mexPrintf("z_loop = %d\n", lz + size_x);

		if (abs(z_diff) < 1e-8f) {
			//if (lo > 3629023)
			//	mexPrintf("lo5 = %d\n", lo);

			//mexPrintf("erotus = %f\n", abs(y_diff));
			//mexPrintf("lo1 = %d\n", lo);

			if (abs(y_diff) < 1e-8) {


				if (yd <= maxyy && yd >= minyy) {
					double minvalue = maxyy * 100.;
					int apu;
					for (int ii = 0; ii < Ny; ii++) {
						double temp = abs(yy_vec[ii] - yd);
						if (temp < minvalue) {
							minvalue = temp;
							apu = ii;
						}
					}
					double templ_ijk = d;
					vector<int> tempk(Nx, apu * Ny + z_loop * Nx * Ny);
					double temp = d * static_cast<double>(Nx);
					temp = 1. / temp;

					if (attenuation_correction == true) {

						double jelppi = 0.;

						for (int iii = 0; iii < tempk.size(); iii++) {
							jelppi += templ_ijk * -atten[tempk[iii] + iii];
						}
						temp = exp(jelppi) * temp;

					}

					float element = static_cast<float>(templ_ijk * temp);

					for (int ii = 0; ii < tempk.size(); ii++) {
						indices.emplace_back((tempk[ii] + ii));
						elements.emplace_back(element);
					}
					lj++;

					lor[lo - 1][0] = lo;
					lor[lo - 1][1] = Nx;
					continue;
				}
			}
			else if (abs(x_diff) < 1e-8f) {

				//mexPrintf("lo = %d\n", lo);
				if (xd <= maxxx && xd >= minxx) {
					double minvalue = maxxx * 100.;
					int apu;
					for (int ii = 0; ii < Nx; ii++) {
						double temp = abs(xx_vec[ii] - xd);
						if (temp < minvalue) {
							minvalue = temp;
							apu = ii;
						}
					}
					double templ_ijk = d;
					vector<int> tempk(Ny, apu + z_loop * Nx * Ny);
					double temp = d * static_cast<double>(Ny);
					temp = 1. / temp;

					if (attenuation_correction == true) {

						double jelppi = 0.;

						for (int iii = 0; iii < tempk.size(); iii++) {
							jelppi += templ_ijk * -atten[tempk[iii] + iii * Ny];
						}
						temp = exp(jelppi) * temp;

					}

					float element = static_cast<float>(templ_ijk * temp);

					for (int ii = 0; ii < tempk.size(); ii++) {
						indices.emplace_back((tempk[ii] + ii * Ny));
						elements.emplace_back(element);
					}
					lj++;

					lor[lo - 1][0] = lo;
					lor[lo - 1][1] = Ny;
					//mexPrintf("lor = %d\n", lor[lo][0]);
					continue;
				}

				continue;
			}

			//if (lo > 3629023)
			//	mexPrintf("lo6 = %d\n", lo);
			//vector<double> tx(Nx + 1, 0);
			//vector<double> ty(Ny + 1, 0);
			//for (int ii = 0; ii < (Nx + 1); ii++) {
			double tx0 = (bx + iij_vec.front() * d - xs) / (x_diff);
			double ty0 = (by + jjk_vec.front() * d - ys) / (y_diff);
			double txe = (bx + iij_vec.back() * d - xs) / (x_diff);
			double tye = (by + jjk_vec.back() * d - ys) / (y_diff);

			//if (lo > 3629023)
			//	mexPrintf("lo7 = %d\n", lo);
				//mexPrintf("tmin = %f\n", ty[ii]);
			//}

			double txmin = min(tx0, txe);
			double txmax = max(tx0, txe);
			double tymin = min(ty0, tye);
			double tymax = max(ty0, tye);

			//if (lo > 3629023)
			//	mexPrintf("lo8 = %d\n", lo);

			double tmin = max(txmin, tymin);
			double tmax = min(txmax, tymax);

			//lor[lo - 1][0] = lo;
			//lor[lo - 1][1] = tmin;

			//if (lo > 3629023)
			//	mexPrintf("lo9 = %d\n", lo);

			//if (lo == 3629023)
			//	mexPrintf("lo = %d\n", lo);

			if (tmin >= tmax)
				continue;

			int imin, imax, jmin, jmax;
			double pxt, pyt;

			//if (lo > 3629023)
			//	mexPrintf("lo10 = %d\n", lo);

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
				//vector<double> tx_n(tx.begin() + imin, tx.begin() + imax);
			}

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
			}

			//if (lo == 67489) {
			//	mexPrintf("jmin = %d\n", jmin);
			//	mexPrintf("jmax = %d\n", jmax);
			//	mexPrintf("imin = %d\n", imin);
			//	mexPrintf("imax = %d\n", imax);
			//}
			
			//mexPrintf("d = %f\n", y_diff);

			//if (lo == 3629023)
			//	mexPrintf("lo2 = %d\n", lo);
			//if (lo > 3629023)
			//	mexPrintf("lo11 = %d\n", lo);

			vector<double> tx_n(imax - imin + 1, 0.);
			vector<double> ty_n(jmax - jmin + 1, 0.);

			if (imax > imin || imin != 0) {
				int tl = 0;
				for (int ii = imin; ii <= imax; ii++) {
					tx_n[tl] = (bx + iij_vec[ii] * d - xs) / (x_diff);
					tl++;
				}
			}
			else
				tx_n[0] = tx0;

			if (jmax > jmin || jmin != 0) {
				int tl = 0;
				for (int ii = jmin; ii <= jmax; ii++) {
					ty_n[tl] = (by + jjk_vec[ii] * d - ys) / (y_diff);
					
					tl++;
				}
			}
			else
				ty_n[0]=ty0;

			//if (lo == 3629023)
			//	mexPrintf("lo3 = %d\n", lo);

			//if (lo > 3629023)
			//	mexPrintf("lo12 = %d\n", lo);

			double L = sqrt(x_diff * x_diff + y_diff * y_diff);

			vector<double> tt;
			tt.reserve(tx_n.size() + ty_n.size() + 1);
			tt.emplace_back(tmin);
			tt.insert(tt.end(), tx_n.begin(), tx_n.end());
			tt.insert(tt.end(), ty_n.begin(), ty_n.end());

			sort(tt.begin(), tt.end());
			tt.erase(unique(tt.begin(), tt.end()), tt.end());

			//if (lo > 3629023)
			//	mexPrintf("lo13 = %d\n", lo);

			//if (lo == 14575) {
			//	for (int kk = 0; kk < tt.size(); kk++)
			//		mexPrintf("tt = %f\n", tt[kk]);
			//	mexPrintf("tt = %d\n", tt.size());
			//}

			int koko = tt.size() - 2;

			int tempi, tempj, tempk;
			//vector<Pair> temp_all;
			//temp_all.reserve(koko + 1);
			vector<double> templ_ijk(koko + 1, 0.);
			vector<int> temp_koko(koko + 1, 0);

			//indices.reserve(tt.size() - 1);
			//elements.reserve(tt.size() - 1);

			tempk = z_loop * Ny * Nx;

			//if (lo > 3629023)
			//	mexPrintf("lo14 = %d\n", lo);

			//if (lo == 67489) {
			//	mexPrintf("tempk = %d\n", tempk);
			//	mexPrintf("z_loop = %d\n", z_loop);
			//}

			double jelppi;

			double temp = 0.;

			//if (lo == 3629023)
			//	mexPrintf("lo4 = %d\n", lo);

			for (int ii = 0; ii <= koko; ii++) {
				jelppi = tt[ii + 1] + tt[ii];
				pxt = xs + (jelppi / 2.) * x_diff;
				tempi = static_cast<int>(floor((pxt - bx) / d));

				pyt = ys + (jelppi / 2.) * y_diff;
				tempj = static_cast<int>(floor((pyt - by) / d));

				temp_koko[ii] = tempj * Ny + tempi + tempk;
				//temp_all[ii].temp_koko.emplace_back(tempj * Ny + tempi + tempk);

				jelppi = tt[ii + 1] - tt[ii];

				templ_ijk[ii] = jelppi * L;
				//temp_all[ii].templ_ijk.emplace_back(jelppi * L);

				temp += templ_ijk[ii];
				//temp += temp_all.at(ii).templ_ijk.front();

				//if (la == 0) {
				//	//mexPrintf("jelppi = %f\n", tt[ii]);
				//	mexPrintf("templ = %f\n", templ_ijk[ii]);
				//	//mexPrintf("L = %f\n", tt[ii + 1]);
				//}
			}

			//if (lo > 3629023)
			//	mexPrintf("lo15 = %d\n", lo);

			temp = 1 / temp;

			//std::sort(make_sort_permute_iter(&temp_koko.begin(), &templ_ijk.begin()),
				//make_sort_permute_iter(&temp_koko.end(), &templ_ijk.end()),	sort_permute_iter_compare());

			//typedef boost::tuple<int&, double&> tup_t;

			//ranges::sort(ranges::view::zip(temp_koko, temp_koko), std::less<>{}, &std::pair<int, int>::first);

			//boost::sort(zip(temp_koko, templ_ijk), [](tup_t i, tup_t j) { return i.get<0>() < j.get<0>(); });

			//std::sort(boost::make_zip_iterator(boost::make_tuple(templ_ijk.begin(), temp_koko.begin())),
			//	boost::make_zip_iterator(boost::make_tuple(templ_ijk.end(), temp_koko.end())));

			//sort(templ_ijk.begin(), templ_ijk.end(), MyComparator(temp_koko));

			//sort(temp_all.begin(), temp_all.end());

			//int minimi = distance(temp_koko.begin(), min_element(temp_koko.begin(), temp_koko.end()));

			//int i = minimi;

			auto i = sort_indexes(temp_koko);

			if (attenuation_correction == true) {

				jelppi = 0.;

				//if (lo == 3629023)
				//	mexPrintf("lo5 = %d\n", lo);

				//if (lo == 3629024)
				//	mexPrintf("koko = %d\n", koko);

				for (int iii = 0; iii <= koko; iii++) {
					jelppi += templ_ijk[iii] * -atten[temp_koko[iii]];
					//jelppi += temp_all.at(iii).templ_ijk.front() * -atten[temp_all.at(iii).temp_koko.front()];
					//if (lo == 3629024) {
					//	mexPrintf("temp_koko = %d\n", atten[temp_koko[iii]]);
					//mexPrintf("i = %d\n", i[iii]);
					//}
				}
				temp = exp(jelppi) * temp;
			}

			for (int ii = 0; ii <= koko; ii++) {


				//if (lo == 3629024)
				//	mexPrintf("ii = %d\n", ii);

				//indices.emplace_back(temp_all.at(ii).temp_koko.front());
				indices.emplace_back(temp_koko[i[ii]]);

				

				//elements.emplace_back(static_cast<float>(temp_all.at(ii).templ_ijk.front() * jelppi));
				elements.emplace_back(static_cast<float>(templ_ijk[i[ii]] * temp));

				

				//if (lo == 34402) {
				//	mexPrintf("koko = %d\n", temp_koko[ii]);
				//	mexPrintf("koko = %d\n", minimi3);
				//	mexPrintf("i = %d\n", i);
				//	//mexPrintf("joku = %d\n", distance(temp_koko.begin(), min_element(temp_koko.begin() + minimi + 1, temp_koko.end())));
				//}

				//if (koko > 1) {

				//	if (ii == 0 && koko > 1 && i < koko && i > 0 && temp_koko[i - 1] > temp_koko[i + 1])
				//		i++;
				//	else if (ii == 0 && koko > 1 && i > 0 && i < koko && temp_koko[i - 1] < temp_koko[i + 1])
				//		i--;
				//	else if (ii == 0 && i == koko && temp_koko[i - 1] < temp_koko[0])
				//		i--;
				//	else if (ii == 0 && i == koko && temp_koko[i - 1] < temp_koko.end()[-1])
				//		i = 0;
				//	else if (ii == 0 && i == koko) {
				//		i = distance(temp_koko.begin(), min_element(temp_koko.begin(), temp_koko.begin() + minimi - 1));
				//		minimi2 = i;
				//	}
				//	else if (ii == 0 && i == 0) {
				//		i = distance(temp_koko.begin(), min_element(temp_koko.begin() + minimi + 1, temp_koko.end()));
				//		minimi3 = i;
				//	}
				//	else if (ii == 0 && i == 0 && temp_koko[i + 1] > temp_koko.end()[-1])
				//		i = koko;
				//	else if (ii == 0 && i == 0) {
				//		i = distance(temp_koko.begin(), min_element(temp_koko.begin() + 1, temp_koko.end()));
				//		minimi3 = i;
				//	}
				//	else if (ii < koko && i == minimi - 1 && minimi2 > 0 || ii < koko && minimi2 > 0 && i == minimi2 - 1 || ii < koko && minimi2_2 > 0 && i == minimi2_2 - 1) {
				//		//mexPrintf("min2 = %d\n", minimi2);
				//		i = distance(temp_koko.begin(), min_element(temp_koko.begin(), temp_koko.begin() - 1 + minimi2));
				//		minimi2_2 = minimi2;
				//		minimi2 = i;
				//	}
				//	else if (ii < koko && i == minimi + 1 && minimi3 > 0 || ii < koko && minimi3 > 0 && i == minimi3 + 1 || ii < koko && minimi3_2 > 0 && i == minimi3_2 + 1) {
				//		//mexPrintf("min3 = %d\n", minimi3);
				//		i = distance(temp_koko.begin(), min_element(temp_koko.begin() + minimi3 + 1, temp_koko.end()));
				//		minimi3_2 = minimi3;
				//		minimi3 = i;
				//	}
				//	else if (i == koko && ii < koko) {
				//		if (minimi == 0 || minimi3 >= 0)
				//			i--;
				//		else {
				//		//mexPrintf("koko = %d\n", temp_koko[ii]);
				//		//mexPrintf("i = %d\n", i);
				//		//mexPrintf("ii = %d\n", ii);
				//		//mexPrintf("joku = %d\n", minimi);
				//			i = distance(temp_koko.begin(), min_element(temp_koko.begin(), temp_koko.begin() + minimi - 1));
				//			minimi2 = i;
				//		}
				//	}
				//	else if (ii > 0 && i == 0 && ii < koko) {
				//		if (minimi == koko || minimi2 >= 0)
				//			i++;
				//		else {
				//			i = distance(temp_koko.begin(), min_element(temp_koko.begin() + minimi + 1, temp_koko.end()));
				//			minimi3 = i;
				//		}
				//	}
				//	else if (i > 0 && temp_koko[i + 1] < temp_koko[i] && temp_koko[i - 1] > temp_koko[i])
				//		i--;
				//	else if (i > 0 && temp_koko[i + 1] > temp_koko[i] && temp_koko[i - 1] < temp_koko[i + 1])
				//		i--;
				//	else
				//		i++;

				//	//if (ii < koko && i == minimi - 1 && minimi2 >= 0)
				//	//	i = distance(temp_koko.begin(), min_element(temp_koko.begin(), temp_koko.begin() - 1 + minimi2));
				//	//else if (ii < koko && i == minimi - 1 && minimi3 >= 0)
				//	//	i = distance(temp_koko.begin(), min_element(temp_koko.begin() + minimi3 + 1, temp_koko.end()));

				//}
				//else {
				//	if (minimi == 0)
				//		i++;
				//	else
				//		i--;
				//}
				//if (indices.end()[-1] == 1137 && indices.end()[-2] == 1136 && indices.end()[-3] == 1238)
				//	mexPrintf("i = %d\n", lo);
				
			}
			//la++;

			lor[lo - 1][0] = lo;
			lor[lo - 1][1] = tt.size() - 1;

			lj++;
			//mexPrintf("lo = %d\n", lor[lo - 1][0]);

			//if (lo > 3629023)
			//	mexPrintf("lo17 = %d\n", lo);

			//if (lo == 3629023)
			//	mexPrintf("lo6 = %d\n", lo);
			//mexPrintf("lo4 = %d\n", lo);
			continue;
		}
		
		else {
			//mexPrintf("lo = %d\n", lo);
			if (abs(y_diff) < 1e-8f) {
				if (yd <= maxyy && yd >= minyy) {

					//vector<double> tx(Nx + 1, 0);
					//vector<double> tz(Nz + 1, 0);
					//for (int ii = 0; ii < (Nz + 1); ii++) {
					//	if (ii < Nx + 1)
					//		tx[ii] = (bx + iij_vec[ii] - xs) / (x_diff);

					//	tz[ii] = (bz + kkj_vec[ii] - zs) / (z_diff);
					//	//mexPrintf("tmin = %f\n", ty[ii]);
					//}

					double tx0 = (bx + iij_vec.front() * d - xs) / (x_diff);
					double tz0 = (bz + kkj_vec.front() * dz - zs) / (z_diff);
					double txe = (bx + iij_vec.back() * d - xs) / (x_diff);
					double tze = (bz + kkj_vec.back() * dz - zs) / (z_diff);

					double txmin = min(tx0, txe);
					double txmax = max(tx0, txe);
					double tzmin = min(tz0, tze);
					double tzmax = max(tz0, tze);

					double tmin = max(txmin, tzmin);
					double tmax = min(txmax, tzmax);

					if (tmin >= tmax)
						continue;

					int imin, imax, kmin, kmax;
					double pxt, pzt;

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
					}

					vector<double> tx_n(imax - imin + 1, 0.);
					vector<double> tz_n(kmax - kmin + 1, 0.);

					if (imax > imin || imin != 0) {
						int tl = 0;
						for (int ii = imin; ii <= imax; ii++) {
							tx_n[tl] = (bx + iij_vec[ii] * d - xs) / (x_diff);
							tl++;
						}
					}
					else
						tx_n[0] = tx0;

					if (kmax > kmin || kmin != 0) {
						int tl = 0;
						for (int ii = kmin; ii <= kmax; ii++) {
							tz_n[tl] = (bz + kkj_vec[ii] * dz - zs) / (z_diff);
							tl++;
						}
					}
					else
						tz_n[0] = tz0;


					double L = sqrt(x_diff * x_diff + z_diff * z_diff);

					vector<double> tt;
					tt.reserve(tx_n.size() + tz_n.size() + 1);
					tt.emplace_back(tmin);
					tt.insert(tt.end(), tx_n.begin(), tx_n.end());
					tt.insert(tt.end(), tz_n.begin(), tz_n.end());

					sort(tt.begin(), tt.end());
					tt.erase(unique(tt.begin(), tt.end()), tt.end());


					int tempi, tempj, tempk;
					double apu2, apu1;

					for (int ii = 0; ii < Ny; ii++) {
						apu1 = abs(yy_vec[ii] - yd);
						if (ii > 0 && apu1 < apu2 || ii == 0) {
							tempj = ii;
							apu2 = apu1;
						}

					}

					tempj = tempj * Ny;

					int koko = tt.size() - 2;

					vector<double> templ_ijk(koko + 1, 0.);
					vector<int> temp_koko(koko + 1, 0);

					double jelppi;

					double temp = 0.;

					for (int ii = 0; ii < tt.size() - 1; ii++) {
						jelppi = (tt[ii + 1] + tt[ii]);

						pxt = xs + (jelppi / 2) * x_diff;
						tempi = static_cast<int>(floor((pxt - bx) / d));

						pzt = zs + (jelppi / 2) * z_diff;
						tempk = static_cast<int>(floor((pzt - bz) / dz));

						temp_koko[ii] = (tempj + tempi + tempk * Ny * Nx);
						//temp_all[ii].temp_koko.emplace_back(tempj * Ny + tempi + tempk);

						jelppi = tt[ii + 1] - tt[ii];

						templ_ijk[ii] = jelppi * L;

						temp += templ_ijk[ii];
					}

					
					//for_each(templ_ijk.begin(), templ_ijk.end(), [&](int n) {
					//	temp += n;
					//});
					temp = 1 / temp;

					//for (int ii = 0; ii < tt.size() - 1; ii++) {
					//	elements.emplace_back(templ_ijk[ii] * temp);
					//	//if (ii == 0) {
					//	//	elements[lo][ii] = templ_ijk[ii] * temp;
					//	//}
					//	//else {
					//	//	elements[lo].emplace_back(templ_ijk[ii] * temp);
					//	//}
					//}

					auto i = sort_indexes(temp_koko);

					//int minimi = distance(temp_koko.begin(), min_element(temp_koko.begin(), temp_koko.end()));

					//int i = minimi;

					//int minimi2 = -1;
					//int minimi2_2 = -1;
					//int minimi3 = -1;
					//int minimi3_2 = -1;

					if (attenuation_correction == true) {

						jelppi = 0.;

						for (int iii = 0; iii <= koko; iii++) {
							jelppi += templ_ijk[iii] * -atten[temp_koko[iii]];
							//jelppi += temp_all.at(iii).templ_ijk.front() * -atten[temp_all.at(iii).temp_koko.front()];
						}
						temp = exp(jelppi) * temp;
					}

					for (int ii = 0; ii <= koko; ii++) {

						//indices.emplace_back(temp_all.at(ii).temp_koko.front());
						indices.emplace_back(temp_koko[i[ii]]);



						//elements.emplace_back(static_cast<float>(temp_all.at(ii).templ_ijk.front() * jelppi));
						elements.emplace_back(static_cast<float>(templ_ijk[i[ii]] * temp));

						//if (lo == 34402) {
						//	mexPrintf("koko = %d\n", temp_koko[ii]);
						//	mexPrintf("koko = %d\n", minimi3);
						//	mexPrintf("i = %d\n", i);
						//	//mexPrintf("joku = %d\n", distance(temp_koko.begin(), min_element(temp_koko.begin() + minimi + 1, temp_koko.end())));
						//}

						//if (koko > 1) {

						//	if (ii == 0 && koko > 1 && i < koko && i > 0 && temp_koko[i - 1] > temp_koko[i + 1])
						//		i++;
						//	else if (ii == 0 && koko > 1 && i > 0 && i < koko && temp_koko[i - 1] < temp_koko[i + 1])
						//		i--;
						//	else if (ii == 0 && i == koko && temp_koko[i - 1] < temp_koko[0])
						//		i--;
						//	else if (ii == 0 && i == koko && temp_koko[i - 1] < temp_koko.end()[-1])
						//		i = 0;
						//	else if (ii == 0 && i == koko) {
						//		i = distance(temp_koko.begin(), min_element(temp_koko.begin(), temp_koko.begin() + minimi - 1));
						//		minimi2 = i;
						//	}
						//	else if (ii == 0 && i == 0) {
						//		i = distance(temp_koko.begin(), min_element(temp_koko.begin() + minimi + 1, temp_koko.end()));
						//		minimi3 = i;
						//	}
						//	else if (ii == 0 && i == 0 && temp_koko[i + 1] > temp_koko.end()[-1])
						//		i = koko;
						//	else if (ii == 0 && i == 0) {
						//		i = distance(temp_koko.begin(), min_element(temp_koko.begin() + 1, temp_koko.end()));
						//		minimi3 = i;
						//	}
						//	else if (ii < koko && i == minimi - 1 && minimi2 > 0 || ii < koko && minimi2 > 0 && i == minimi2 - 1 || ii < koko && minimi2_2 > 0 && i == minimi2_2 - 1) {
						//		//mexPrintf("min2 = %d\n", minimi2);
						//		i = distance(temp_koko.begin(), min_element(temp_koko.begin(), temp_koko.begin() - 1 + minimi2));
						//		minimi2_2 = minimi2;
						//		minimi2 = i;
						//	}
						//	else if (ii < koko && i == minimi + 1 && minimi3 > 0 || ii < koko && minimi3 > 0 && i == minimi3 + 1 || ii < koko && minimi3_2 > 0 && i == minimi3_2 + 1) {
						//		//mexPrintf("min3 = %d\n", minimi3);
						//		i = distance(temp_koko.begin(), min_element(temp_koko.begin() + minimi3 + 1, temp_koko.end()));
						//		minimi3_2 = minimi3;
						//		minimi3 = i;
						//	}
						//	else if (i == koko && ii < koko) {
						//		if (minimi == 0 || minimi3 >= 0)
						//			i--;
						//		else {
						//			//mexPrintf("koko = %d\n", temp_koko[ii]);
						//			//mexPrintf("i = %d\n", i);
						//			//mexPrintf("ii = %d\n", ii);
						//			//mexPrintf("joku = %d\n", minimi);
						//			i = distance(temp_koko.begin(), min_element(temp_koko.begin(), temp_koko.begin() + minimi - 1));
						//			minimi2 = i;
						//		}
						//	}
						//	else if (ii > 0 && i == 0 && ii < koko) {
						//		if (minimi == koko || minimi2 >= 0)
						//			i++;
						//		else {
						//			i = distance(temp_koko.begin(), min_element(temp_koko.begin() + minimi + 1, temp_koko.end()));
						//			minimi3 = i;
						//		}
						//	}
						//	else if (i > 0 && temp_koko[i + 1] < temp_koko[i] && temp_koko[i - 1] > temp_koko[i])
						//		i--;
						//	else if (i > 0 && temp_koko[i + 1] > temp_koko[i] && temp_koko[i - 1] < temp_koko[i + 1])
						//		i--;
						//	else
						//		i++;

						//	//if (ii < koko && i == minimi - 1 && minimi2 >= 0)
						//	//	i = distance(temp_koko.begin(), min_element(temp_koko.begin(), temp_koko.begin() - 1 + minimi2));
						//	//else if (ii < koko && i == minimi - 1 && minimi3 >= 0)
						//	//	i = distance(temp_koko.begin(), min_element(temp_koko.begin() + minimi3 + 1, temp_koko.end()));

						//}
						//else {
						//	if (minimi == 0)
						//		i++;
						//	else
						//		i--;
						//}
						//if (indices.end()[-1] == 1137 && indices.end()[-2] == 1136 && indices.end()[-3] == 1238)
						//	mexPrintf("i = %d\n", lo);

					}


					lor[lo - 1][0] = lo;
					lor[lo - 1][1] = tt.size() - 1;
					lj++;
					continue;

				}
				continue;
			}
			else if (abs(x_diff) < 1e-8f) {
				if (xd <= maxxx && xd >= minxx) {

					double ty0 = (by + jjk_vec.front() * d - ys) / (y_diff);
					double tz0 = (bz + kkj_vec.front() * dz - zs) / (z_diff);
					double tye = (by + jjk_vec.back() * d - ys) / (y_diff);
					double tze = (bz + kkj_vec.back() * dz - zs) / (z_diff);

					double tymin = min(ty0, tye);
					double tymax = max(ty0, tye);
					double tzmin = min(tz0, tze);
					double tzmax = max(tz0, tze);

					double tmin = max(tymin, tzmin);
					double tmax = min(tymax, tzmax);

					if (tmin >= tmax)
						continue;


					int jmin, jmax, kmin, kmax;
					double pyt, pzt;

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
					}

					vector<double> ty_n(jmax - jmin + 1, 0.);
					vector<double> tz_n(kmax - kmin + 1, 0.);

					if (jmax > jmin || jmin != 0) {
						int tl = 0;
						for (int ii = jmin; ii <= jmax; ii++) {
							ty_n[tl] = (by + jjk_vec[ii] * d - ys) / (y_diff);
							tl++;
						}
					}
					else
						ty_n[0] = ty0;

					if (kmax > kmin || kmin != 0) {
						int tl = 0;
						for (int ii = kmin; ii <= kmax; ii++) {
							tz_n[tl] = (bz + kkj_vec[ii] * dz - zs) / (z_diff);
							tl++;
						}
					}
					else
						tz_n[0] = tz0;


					double L = sqrt(y_diff * y_diff + z_diff * z_diff);

					vector<double> tt;
					tt.reserve(ty_n.size() + tz_n.size() + 1);
					tt.emplace_back(tmin);
					tt.insert(tt.end(), ty_n.begin(), ty_n.end());
					tt.insert(tt.end(), tz_n.begin(), tz_n.end());

					sort(tt.begin(), tt.end());
					tt.erase(unique(tt.begin(), tt.end()), tt.end());

					int tempi, tempj, tempk;
					double apu2, apu1;

					for (int ii = 0; ii < Nx; ii++) {
						apu1 = abs(xx_vec[ii] - xd);
						if (ii > 0 && apu1 < apu2 || ii == 0) {
							tempi = ii;
							apu2 = apu1;
						}

					}

					//tempi = tempi * Ny;

					int koko = tt.size() - 2;

					vector<double> templ_ijk(koko + 1, 0.);
					vector<int> temp_koko(koko + 1, 0);

					double jelppi;

					double temp = 0.;

					for (int ii = 0; ii < tt.size() - 1; ii++) {
						jelppi = (tt[ii + 1] + tt[ii]);

						pyt = ys + (jelppi / 2) * y_diff;
						tempj = static_cast<int>(floor((pyt - by) / d));

						pzt = zs + (jelppi / 2) * z_diff;
						tempk = static_cast<int>(floor((pzt - bz) / dz));

						temp_koko[ii] = (tempj * Ny + tempi + tempk * Ny * Nx);

						jelppi = tt[ii + 1] - tt[ii];

						templ_ijk[ii] = jelppi * L;

						temp += templ_ijk[ii];

						//mexPrintf("ii = %d\n", (tempj * Ny + tempi + tempk * Ny * Nx));

						//if (ii == 0) {
						//	indices[lo][ii] = (int)(tempj * Ny + tempi + tempk * Ny * Nx);
						//}
						//else {
						//	indices[lo].emplace_back((int)(tempj * Ny + tempi + tempk * Ny * Nx));
						//}

					}

					
					//for_each(templ_ijk.begin(), templ_ijk.end(), [&](int n) {
					//	temp += n;
					//});
					temp = 1 / temp;

					auto i = sort_indexes(temp_koko);

					//int minimi = distance(temp_koko.begin(), min_element(temp_koko.begin(), temp_koko.end()));

					//int i = minimi;

					//int minimi2 = -1;
					//int minimi2_2 = -1;
					//int minimi3 = -1;
					//int minimi3_2 = -1;

					if (attenuation_correction == true) {

						jelppi = 0.;

						for (int iii = 0; iii <= koko; iii++) {
							jelppi += templ_ijk[iii] * -atten[temp_koko[iii]];
							//jelppi += temp_all.at(iii).templ_ijk.front() * -atten[temp_all.at(iii).temp_koko.front()];
						}
						temp = exp(jelppi) * temp;

					}

					for (int ii = 0; ii <= koko; ii++) {

						

						//indices.emplace_back(temp_all.at(ii).temp_koko.front());
						indices.emplace_back(temp_koko[i[ii]]);



						//elements.emplace_back(static_cast<float>(temp_all.at(ii).templ_ijk.front() * jelppi));
						elements.emplace_back(static_cast<float>(templ_ijk[i[ii]] * temp));

						//if (lo == 34402) {
						//	mexPrintf("koko = %d\n", temp_koko[ii]);
						//	mexPrintf("koko = %d\n", minimi3);
						//	mexPrintf("i = %d\n", i);
						//	//mexPrintf("joku = %d\n", distance(temp_koko.begin(), min_element(temp_koko.begin() + minimi + 1, temp_koko.end())));
						//}

						//if (koko > 1) {

						//	if (ii == 0 && koko > 1 && i < koko && i > 0 && temp_koko[i - 1] > temp_koko[i + 1])
						//		i++;
						//	else if (ii == 0 && koko > 1 && i > 0 && i < koko && temp_koko[i - 1] < temp_koko[i + 1])
						//		i--;
						//	else if (ii == 0 && i == koko && temp_koko[i - 1] < temp_koko[0])
						//		i--;
						//	else if (ii == 0 && i == koko && temp_koko[i - 1] < temp_koko.end()[-1])
						//		i = 0;
						//	else if (ii == 0 && i == koko) {
						//		i = distance(temp_koko.begin(), min_element(temp_koko.begin(), temp_koko.begin() + minimi - 1));
						//		minimi2 = i;
						//	}
						//	else if (ii == 0 && i == 0) {
						//		i = distance(temp_koko.begin(), min_element(temp_koko.begin() + minimi + 1, temp_koko.end()));
						//		minimi3 = i;
						//	}
						//	else if (ii == 0 && i == 0 && temp_koko[i + 1] > temp_koko.end()[-1])
						//		i = koko;
						//	else if (ii == 0 && i == 0) {
						//		i = distance(temp_koko.begin(), min_element(temp_koko.begin() + 1, temp_koko.end()));
						//		minimi3 = i;
						//	}
						//	else if (ii < koko && i == minimi - 1 && minimi2 > 0 || ii < koko && minimi2 > 0 && i == minimi2 - 1 || ii < koko && minimi2_2 > 0 && i == minimi2_2 - 1) {
						//		//mexPrintf("min2 = %d\n", minimi2);
						//		i = distance(temp_koko.begin(), min_element(temp_koko.begin(), temp_koko.begin() - 1 + minimi2));
						//		minimi2_2 = minimi2;
						//		minimi2 = i;
						//	}
						//	else if (ii < koko && i == minimi + 1 && minimi3 > 0 || ii < koko && minimi3 > 0 && i == minimi3 + 1 || ii < koko && minimi3_2 > 0 && i == minimi3_2 + 1) {
						//		//mexPrintf("min3 = %d\n", minimi3);
						//		i = distance(temp_koko.begin(), min_element(temp_koko.begin() + minimi3 + 1, temp_koko.end()));
						//		minimi3_2 = minimi3;
						//		minimi3 = i;
						//	}
						//	else if (i == koko && ii < koko) {
						//		if (minimi == 0 || minimi3 >= 0)
						//			i--;
						//		else {
						//			//mexPrintf("koko = %d\n", temp_koko[ii]);
						//			//mexPrintf("i = %d\n", i);
						//			//mexPrintf("ii = %d\n", ii);
						//			//mexPrintf("joku = %d\n", minimi);
						//			i = distance(temp_koko.begin(), min_element(temp_koko.begin(), temp_koko.begin() + minimi - 1));
						//			minimi2 = i;
						//		}
						//	}
						//	else if (ii > 0 && i == 0 && ii < koko) {
						//		if (minimi == koko || minimi2 >= 0)
						//			i++;
						//		else {
						//			i = distance(temp_koko.begin(), min_element(temp_koko.begin() + minimi + 1, temp_koko.end()));
						//			minimi3 = i;
						//		}
						//	}
						//	else if (i > 0 && temp_koko[i + 1] < temp_koko[i] && temp_koko[i - 1] > temp_koko[i])
						//		i--;
						//	else if (i > 0 && temp_koko[i + 1] > temp_koko[i] && temp_koko[i - 1] < temp_koko[i + 1])
						//		i--;
						//	else
						//		i++;

						//	//if (ii < koko && i == minimi - 1 && minimi2 >= 0)
						//	//	i = distance(temp_koko.begin(), min_element(temp_koko.begin(), temp_koko.begin() - 1 + minimi2));
						//	//else if (ii < koko && i == minimi - 1 && minimi3 >= 0)
						//	//	i = distance(temp_koko.begin(), min_element(temp_koko.begin() + minimi3 + 1, temp_koko.end()));

						//}
						//else {
						//	if (minimi == 0)
						//		i++;
						//	else
						//		i--;
						//}
						//if (indices.end()[-1] == 1137 && indices.end()[-2] == 1136 && indices.end()[-3] == 1238)
						//	mexPrintf("i = %d\n", lo);

					}


					lor[lo - 1][0] = lo;
					lor[lo - 1][1] = tt.size() - 1;
					lj++;
					continue;

				}
				continue;
			}

			

			double ty0 = (by + jjk_vec.front() * d - ys) / (y_diff);
			double tz0 = (bz + kkj_vec.front() * dz - zs) / (z_diff);
			double tye = (by + jjk_vec.back() * d - ys) / (y_diff);
			double tze = (bz + kkj_vec.back() * dz - zs) / (z_diff);
			double tx0 = (bx + iij_vec.front() * d - xs) / (x_diff);
			double txe = (bx + iij_vec.back() * d - xs) / (x_diff);

			double txmin = min(tx0, txe);
			double txmax = max(tx0, txe);
			double tymin = min(ty0, tye);
			double tymax = max(ty0, tye);
			double tzmin = min(tz0, tze);
			double tzmax = max(tz0, tze);

			double tmin = max(max(txmin, tzmin), tymin);
			double tmax = min(min(txmax, tzmax), tymax);

			//if (lo == 13877283) {
			//	mexPrintf("tmin = %f\n", tmin);
			//	mexPrintf("tmax = %f\n", tmax);
			//}

			//if (lo == 13877283)
			//	mexPrintf("lo1 = %d\n", lo);

			if (tmin >= tmax)
				continue;

			int imin, imax, jmin, jmax, kmin, kmax;
			double pxt, pyt, pzt;


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
			}

			//if (lo == 13103601) {
			//	mexPrintf("kmin = %f\n", kmin);
			//	mexPrintf("kmax = %f\n", kmax);
			//}

			//if (lo == 13877283)
			//	mexPrintf("lo2 = %d\n", lo);

			vector<double> tx_n(imax - imin + 1, 0.);
			vector<double> tz_n(kmax - kmin + 1, 0.);
			vector<double> ty_n(jmax - jmin + 1, 0.);

			if (jmax > jmin || jmin != 0) {
				int tl = 0;
				for (int ii = jmin; ii <= jmax; ii++) {
					ty_n[tl] = (by + jjk_vec[ii] * d - ys) / (y_diff);
					//if (lo == 13103601)
					//	mexPrintf("ty = %f\n", ty_n[tl]);
					tl++;
				}
			}
			else
				ty_n[0] = ty0;

			if (imax > imin || imin != 0) {
				int tl = 0;
				for (int ii = imin; ii <= imax; ii++) {
					tx_n[tl] = (bx + iij_vec[ii] * d - xs) / (x_diff);
					//if (lo == 13103601)
					//	mexPrintf("tx = %f\n", tx_n[tl]);
					tl++;
				}
			}
			else
				tx_n[0] = tx0;

			if (kmax > kmin || kmin != 0) {
				int tl = 0;
				for (int ii = kmin; ii <= kmax; ii++) {
					tz_n[tl] = (bz + kkj_vec[ii] * dz - zs) / (z_diff);
					//if (lo == 13103601)
					//	mexPrintf("tz = %f\n", bz + kkj_vec[ii] * dz);
					tl++;
				}
			}
			else
				tz_n[0] = tz0;

			double L = sqrt(x_diff * x_diff + z_diff * z_diff + y_diff * y_diff);

			//if (lo == 13877283)
			//	mexPrintf("lo3 = %d\n", lo);

			vector<double> tt;
			tt.reserve(tx_n.size() + tz_n.size() + ty_n.size() + 1);
			tt.emplace_back(tmin);
			tt.insert(tt.end(), tx_n.begin(), tx_n.end());
			tt.insert(tt.end(), ty_n.begin(), ty_n.end());
			tt.insert(tt.end(), tz_n.begin(), tz_n.end());

			sort(tt.begin(), tt.end());
			tt.erase(unique(tt.begin(), tt.end()), tt.end());

			//if (lo == 36976102) {
			//	for (int ll = 0; ll < tt.size(); ll++) {
			//		mexPrintf("tt = %f\n", tt[ll]);
			//	}
			//}

			int tempi, tempj, tempk;

			int koko = tt.size() - 2;

			//if (lo == 13877283)
			//	mexPrintf("lo4 = %d\n", lo);

			vector<double> templ_ijk(koko + 1, 0.);
			vector<int> temp_koko(koko + 1, 0);

			//if (lo == 13877283)
			//	mexPrintf("lo5 = %d\n", lo);

			double jelppi;

			double temp = 0.;

			for (int ii = 0; ii < tt.size() - 1; ii++) {
				jelppi = (tt[ii + 1] + tt[ii]);

				pyt = ys + (jelppi / 2) * y_diff;
				tempj = static_cast<int>(floor((pyt - by) / d));

				pzt = zs + (jelppi / 2) * z_diff;
				tempk = static_cast<int>(floor((pzt - bz) / dz));

				pxt = xs + (jelppi / 2) * x_diff;
				tempi = static_cast<int>(floor((pxt - bx) / d));

				temp_koko[ii] = (tempj * Ny + tempi + tempk * Ny * Nx);

				jelppi = tt[ii + 1] - tt[ii];

				templ_ijk[ii] = jelppi * L;

				temp += templ_ijk[ii];

				//if (lo == 36976102)
				//	mexPrintf("templ = %f\n", templ_ijk[ii]);

				//mexPrintf("ii = %d\n", (tempj * Ny + tempi + tempk * Ny * Nx));

				//if (ii == 0) {
				//	indices[lo][ii] = (int)(tempj * Ny + tempi + tempk * Ny * Nx);
				//}
				//else {
				//	indices[lo].emplace_back((int)(tempj * Ny + tempi + tempk * Ny * Nx));
				//}

			}

			//if (lo == 13877283)
			//	mexPrintf("lo6 = %d\n", lo);

			//double temp;
			//for_each(templ_ijk.begin(), templ_ijk.end(), [&](int n) {
			//	temp += n;
			//});
			temp = 1 / temp;

			

			//int minimi = distance(temp_koko.begin(), min_element(temp_koko.begin(), temp_koko.end()));

			//if (lo == 36976102)
				//mexPrintf("temp = %f\n", temp);

			//int i = minimi;

			//int minimi2 = -1;
			//int minimi2_2 = -1;
			//int minimi3 = -1;
			//int minimi3_2 = -1;

			if (attenuation_correction == true) {

				jelppi = 0.;

				for (int iii = 0; iii <= koko; iii++) {
					jelppi = jelppi + templ_ijk[iii] * -atten[temp_koko[iii]];
					//if (lo == 36976102) {
					//	mexPrintf("atten = %f\n", jelppi);
					//	mexPrintf("atten_vali = %f\n", templ_ijk[iii] * -atten[temp_koko[iii]]);
					//	mexPrintf("templ = %f\n", templ_ijk[iii]);
					//	mexPrintf("temp_koko = %d\n", temp_koko[iii]);
					//}
					//jelppi += temp_all.at(iii).templ_ijk.front() * -atten[temp_all.at(iii).temp_koko.front()];
				}
				//if (lo == 36976102)
				//mexPrintf("jelppi1 = %f\n", jelppi);
				temp = exp(jelppi) * temp;

			}
			
			auto i = sort_indexes(temp_koko);

			//if (lo == 13877283)
			//	continue;
				//mexPrintf("koko = %d\n", koko);
				//mexPrintf("lo8 = %d\n", lo);

			for (int ii = 0; ii <= koko; ii++) {

				//if (ii == 0) {
					
					//if (lo == 36976102)
						//mexPrintf("jelppi2 = %f\n", jelppi);
					//jelppi = temp;
					//if (la == 0) {
					//	mexPrintf("jelppi = %f\n", jelppi);
					//	mexPrintf("temp = %f\n", temp);
					//	//mexPrintf("L = %f\n", tt[ii + 1]);
					//}
				//}

				//indices.emplace_back(temp_all.at(ii).temp_koko.front());
				indices.emplace_back(temp_koko[i[ii]]);



				//elements.emplace_back(static_cast<float>(temp_all.at(ii).templ_ijk.front() * jelppi));
				elements.emplace_back(static_cast<float>(templ_ijk[i[ii]] * temp));

				//if (lo == 36976102) {
				//	mexPrintf("temp_koko = %f\n", (templ_ijk[i[ii]] * jelppi));
					
				//	mexPrintf("minimi2 = %d\n", minimi2);
				//	mexPrintf("i = %d\n", i);
				//	//mexPrintf("joku = %d\n", distance(temp_koko.begin(), min_element(temp_koko.begin() + minimi + 1, temp_koko.end())));
				//}

				//if (koko > 1) {

				//	if (ii == 0 && koko > 1 && i < koko && i > 0 && temp_koko[i - 1] > temp_koko[i + 1])
				//		i++;
				//	else if (ii == 0 && koko > 1 && i > 0 && i < koko && temp_koko[i - 1] < temp_koko[i + 1])
				//		i--;
				//	else if (ii == 0 && i == koko && temp_koko[i - 1] < temp_koko[0])
				//		i--;
				//	else if (ii == 0 && i == koko && temp_koko[i - 1] < temp_koko.end()[-1])
				//		i = 0;
				//	else if (ii == 0 && i == koko) {
				//		i = distance(temp_koko.begin(), min_element(temp_koko.begin(), temp_koko.begin() + minimi - 1));
				//		minimi2 = i;
				//	}
				//	else if (ii == 0 && i == 0) {
				//		i = distance(temp_koko.begin(), min_element(temp_koko.begin() + minimi + 1, temp_koko.end()));
				//		minimi3 = i;
				//	}
				//	else if (ii == 0 && i == 0 && temp_koko[i + 1] > temp_koko.end()[-1])
				//		i = koko;
				//	else if (ii == 0 && i == 0) {
				//		i = distance(temp_koko.begin(), min_element(temp_koko.begin() + 1, temp_koko.end()));
				//		minimi3 = i;
				//	}
				//	else if (ii < koko && i == minimi - 1 && minimi2 > 0 || ii < koko && minimi2 > 0 && i == minimi2 - 1 || ii < koko && minimi2_2 > 0 && i == minimi2_2 - 1) {
				//		//mexPrintf("min2 = %d\n", minimi2);
				//		i = distance(temp_koko.begin(), min_element(temp_koko.begin(), temp_koko.begin() - 1 + minimi2));
				//		minimi2_2 = minimi2;
				//		minimi2 = i;
				//	}
				//	else if (ii < koko && i == minimi + 1 && minimi3 > 0 || ii < koko && minimi3 > 0 && i == minimi3 + 1 || ii < koko && minimi3_2 > 0 && i == minimi3_2 + 1) {
				//		//mexPrintf("min3 = %d\n", minimi3);
				//		i = distance(temp_koko.begin(), min_element(temp_koko.begin() + minimi3 + 1, temp_koko.end()));
				//		minimi3_2 = minimi3;
				//		minimi3 = i;
				//	}
				//	else if (i == koko && ii < koko) {
				//		if (minimi == 0 || minimi3 >= 0)
				//			i--;
				//		else {
				//			//mexPrintf("koko = %d\n", temp_koko[ii]);
				//			//mexPrintf("i = %d\n", i);
				//			//mexPrintf("ii = %d\n", ii);
				//			//mexPrintf("joku = %d\n", minimi);
				//			i = distance(temp_koko.begin(), min_element(temp_koko.begin(), temp_koko.begin() + minimi - 1));
				//			minimi2 = i;
				//		}
				//	}
				//	else if (ii > 0 && i == 0 && ii < koko) {
				//		if (minimi == koko || minimi2 >= 0)
				//			i++;
				//		else {
				//			i = distance(temp_koko.begin(), min_element(temp_koko.begin() + minimi + 1, temp_koko.end()));
				//			minimi3 = i;
				//		}
				//	}
				//	else if (i > 0 && temp_koko[i + 1] < temp_koko[i] && temp_koko[i - 1] > temp_koko[i])
				//		i--;
				//	else if (i > 0 && temp_koko[i + 1] > temp_koko[i] && temp_koko[i - 1] < temp_koko[i + 1])
				//		i--;
				//	else if (i > 0 && temp_koko[i + 1] < temp_koko[i] && temp_koko[i - 1] > temp_koko[i + 1])
				//		i--;
				//	else
				//		i++;

				//	//if (ii < koko && i == minimi - 1 && minimi2 >= 0)
				//	//	i = distance(temp_koko.begin(), min_element(temp_koko.begin(), temp_koko.begin() - 1 + minimi2));
				//	//else if (ii < koko && i == minimi - 1 && minimi3 >= 0)
				//	//	i = distance(temp_koko.begin(), min_element(temp_koko.begin() + minimi3 + 1, temp_koko.end()));

				//}
				//else {
				//	if (minimi == 0)
				//		i++;
				//	else
				//		i--;
				//}
				//if (indices.end()[-1] == 440640 && indices.end()[-2] == 440640 && indices.end()[-3] == 440704 && indices.end()[-4] == 440640 && indices.end()[-5] == 440704 && indices.end()[-6] == 440640 && indices.end()[-7] == 438309)
				//	mexPrintf("i = %d\n", lo);

			}
			//if (lo == 13877283)
			//	mexPrintf("lo9 = %d\n", lo);

			lor[lo - 1][0] = lo;
			lor[lo - 1][1] = tt.size() - 1;
			lj++;
			continue;

		}
		
	}
	
	//vector<int> global_vector;

	//mexPrintf("lo = %d\n", lor[2][0]);
	//mexPrintf("lo = %d\n", lor[loop_var_par - 1][1]);

	size_t outSize1 = lj * 2;
	size_t outSize2 = 1;

	plhs[0] = mxCreateNumericMatrix(loop_var_par * 2, outSize2, mxINT32_CLASS, mxREAL);

	int* outputMatrix = (int *)mxGetData(plhs[0]);

	plhs[1] = mxCreateNumericMatrix(indices.size(), outSize2, mxINT32_CLASS, mxREAL);

	int* outputMatrix2 = (int *)mxGetData(plhs[1]);

	plhs[2] = mxCreateNumericMatrix(elements.size(), outSize2, mxSINGLE_CLASS, mxREAL);

	float* outputMatrix3 = (float *)mxGetData(plhs[2]);

	//plhs[0] = mxCreateNumericMatrix(outSize1, outSize2, mxSINGLE_CLASS, mxREAL);

	//double* outputMatrix = (double *)mxGetData(plhs[0]);

	//int* outputMatrix = lor.data();

	

	int lb = 0;
	for (int col = 0; col < 2; col++) {
		int la = 0;
		for (int row = 0; row < loop_var_par; row++) {
			outputMatrix[la + col*loop_var_par] = lor[row][col];
			la++;
		}
	}

	lor.erase(lor.begin(), lor.end());
	lor.shrink_to_fit();

	copy(indices.begin(), indices.end(), outputMatrix2);
	indices.erase(indices.begin(), indices.end());
	indices.shrink_to_fit();
	copy(elements.begin(), elements.end(), outputMatrix3);

	//outputMatrix2 = &indices[0];
	//outputMatrix3 = &elements[0];

	//for (auto & buffer : lor) {
	//	move(buffer.begin(), buffer.end(), back_inserter(global_vector));
	//}

	

	//std::memcpy(mxGetData(plhs[0]), outputMatrix, (outSize1)*(outSize2) * sizeof(int));
}