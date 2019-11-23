/**************************************************************************
* This function is used to check for the number of voxels each LOR/ray
* traverses and also whether the LOR/ray actually intersects with the pixel
* space.
* Raw list-mode data and sinogram data have slightly different versions.
* Output is the number of voxels the ray has traversed (if the LOR does not
* traverse the pixel space, this value will be 0).
*
* Copyright (C) 2019  Ville-Veikko Wettenhovi
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
#include "projector_functions.h"

constexpr int TYPE = 0;


using namespace std;

void improved_siddon(const size_t loop_var_par, const uint32_t size_x, const double zmax, const uint32_t TotSinos, uint16_t* lor, const double maxyy,
	const double maxxx, const vector<double>& xx_vec, const vector<double>& z_det_vec, const double dy, const vector<double>& yy_vec,
	const double* x, const double* y, const double* z_det, const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const double dx, const double dz,
	const double bx, const double by, const double bz, const uint32_t block1, const uint32_t blocks, const uint16_t* L, const uint32_t* pseudos,
	const bool raw, const uint32_t pRows, const uint32_t det_per_ring, const uint32_t type, uint16_t* lor_orth, const double crystal_size, const double crystal_size_z, 
	const double* x_center,	const double* y_center, const double* z_center) {


	const double bzb = bz + static_cast<double>(Nz) * dz;

	const int32_T dec = static_cast<int32_T>(ceil(crystal_size_z / sqrt(dz * dz * 2.))) * 3;
	const int32_T decx = static_cast<int32_T>(ceil(crystal_size_z / sqrt(dx * dx * 2.))) * 1;


	ThreadPool::ParallelFor(static_cast<size_t>(0), loop_var_par, [&](uint32_t lo) {
		Det detectors;

		if (raw) {
			const uint32_t detektorit1 = static_cast<uint32_t>(L[lo * 2u]);
			const uint32_t detektorit2 = static_cast<uint32_t>(L[lo * 2u + 1u]);

			if (detektorit1 == detektorit2)
				return;

			const uint32_t loop1 = 1 + ((detektorit1 - 1u) / det_per_ring);
			const uint32_t loop2 = 1 + ((detektorit2 - 1u) / det_per_ring);

			if (loop1 == loop2) {
				if (loop1 > blocks || loop1 < block1 || loop2 > blocks || loop2 < block1)
					return;
				detectors.zs = z_det_vec[loop1 - 1u];
				detectors.zd = detectors.zs;
			}
			else {
				if (loop1 > blocks && loop2 > blocks || loop1 < block1 && loop2 < block1)
					return;
				detectors.zs = z_det_vec[loop1 - 1u];
				detectors.zd = z_det_vec[loop2 - 1u];
			}

			detectors.xs = x[detektorit1 - det_per_ring * (loop1 - 1u) - 1u];
			detectors.xd = x[detektorit2 - det_per_ring * (loop2 - 1u) - 1u];
			detectors.ys = y[detektorit1 - det_per_ring * (loop1 - 1u) - 1u];
			detectors.yd = y[detektorit2 - det_per_ring * (loop2 - 1u) - 1u];
		}
		else {
			const uint32_t id = lo % size_x;
			const uint32_t idz = lo / size_x;

			detectors.xs = x[id];
			detectors.xd = x[id + size_x];
			detectors.ys = y[id];
			detectors.yd = y[id + size_x];
			detectors.zs = z_det[idz];
			detectors.zd = z_det[idz + TotSinos];
		}

		if (detectors.zd < detectors.zs) {
			double tempa = detectors.zd;
			detectors.zd = detectors.zs;
			detectors.zs = tempa;
			tempa = detectors.xd;
			detectors.xd = detectors.xs;
			detectors.xs = tempa;
			tempa = detectors.yd;
			detectors.yd = detectors.ys;
			detectors.ys = tempa;
		}


		const double y_diff = (detectors.yd - detectors.ys);
		const double x_diff = (detectors.xd - detectors.xs);
		const double z_diff = (detectors.zd - detectors.zs);

		double kerroinz = 0., kerroin, length_;

		uint16_t temp_koko_orth = 0u;
		uint16_t temp_koko = 0u;
		uint16_t temp_koko_orth_3D = 0u;

		uint32_t xyz = 0u;
		uint32_t Np = 0u;

		if (type > 0u) {
			if (type == 2u) {
				kerroinz = norm(x_diff, y_diff, z_diff) * crystal_size_z;
			}
			kerroin = detectors.xd * detectors.ys - detectors.yd * detectors.xs;
			length_ = sqrt(y_diff * y_diff + x_diff * x_diff) * crystal_size;
		}

		if (fabs(z_diff) < 1e-8) {

			const uint32_t tempk = z_ring(zmax, detectors.zs, static_cast<double>(NSlices));

			if (fabs(y_diff) < 1e-8) {


				if (detectors.yd <= maxyy && detectors.yd >= by) {
					// The number of voxels the LOR/ray traverses

					if (type > 0u) {
						orth_perpendicular_precompute(Nx, Ny, detectors.yd, yy_vec, -x_diff, y_center, kerroin, length_, temp_koko_orth);
						lor_orth[lo] = temp_koko_orth;
						if (type == 2u) {
							double temppi = detectors.xs;
							detectors.xs = detectors.ys;
							detectors.ys = temppi;
							orth_perpendicular_precompute_3D(Nx, Ny, Nz, detectors.yd, yy_vec, y_center, x_center[0], z_center, kerroinz, temp_koko_orth_3D, 
								detectors, y_diff, x_diff, z_diff, tempk);
							lor_orth[lo + loop_var_par] = temp_koko_orth_3D;
						}
					}
					lor[lo] = static_cast<uint16_t>(Nx);
				}
				// LOR/ray doesn't traverse through the pixel space
			}
			else if (fabs(x_diff) < 1e-8) {

				if (detectors.xd <= maxxx && detectors.xd >= bx) {
					// The number of voxels the LOR/ray traverses
					if (type > 0u) {
						orth_perpendicular_precompute(Ny, Nx, detectors.xd, xx_vec, y_diff, x_center, kerroin, length_, temp_koko_orth);
						lor_orth[lo] = temp_koko_orth;
						if (type == 2u) {
							orth_perpendicular_precompute_3D(Ny, Nx, Nz, detectors.xd, xx_vec, x_center, y_center[0], z_center, kerroinz, temp_koko_orth_3D,
								detectors, x_diff, y_diff, z_diff, tempk);
							lor_orth[lo + loop_var_par] = temp_koko_orth_3D;
						}
					}
					lor[lo] = static_cast<uint16_t>(Ny);
				}
				// LOR/ray doesn't traverse through the pixel space
			}
			else {
				int32_t tempi = 0, tempj = 0, ju = 0, iu = 0;
				double tyu = 0., txu = 0., tc = 0., ty0 = 0., tx0 = 0.;

				const bool skip = siddon_pre_loop_2D(bx, by, x_diff, y_diff, maxxx, maxyy, dx, dy, Nx, Ny, tempi, tempj, txu, tyu, Np, TYPE,
					detectors.ys, detectors.xs, detectors.yd, detectors.xd, tc, iu, ju, tx0, ty0);

				if (!skip) {

					//bool var = false;
					int32_t n_tempk = tempk;
					bool z_bool = true;

					for (uint32_t ii = 0u; ii < Np; ii++) {

						temp_koko++;

						if (tx0 < ty0) {
							if (ii == (Np - 1u) && type > 0u) {
								orth_distance_precompute(tempi, Nx, y_diff, x_diff, x_center, y_center[tempj], kerroin, length_, temp_koko_orth);
								if (type == 2u) {
									orth_distance_precompute_3D(tempi, Nx, Nz, y_diff, x_diff, z_diff, x_center, y_center[tempj], z_center, kerroinz,
										detectors, temp_koko_orth_3D, tempk, lo, n_tempk, dec, decx);
								}
							}
							tempi += iu;
							tx0 += txu;
							xyz = 1u;
						}
						else {
							if (type > 0u) {
								orth_distance_precompute(tempi, Nx, y_diff, x_diff, x_center, y_center[tempj], kerroin, length_, temp_koko_orth);
								if (type == 2u) {
									orth_distance_precompute_3D(tempi, Nx, Nz, y_diff, x_diff, z_diff, x_center, y_center[tempj], z_center, kerroinz,
										detectors, temp_koko_orth_3D, tempk, lo, n_tempk, dec, decx);
								}
							}
							tempj += ju;
							ty0 += tyu;
							xyz = 2u;
						}
						if (tempj < 0 || tempi < 0 || tempi >= Nx || tempj >= Ny) {
							if (xyz == 1u && type > 0u && ii != (Np - 1u)) {
								tempi -= iu;
								orth_distance_precompute(tempi, Nx, y_diff, x_diff, x_center, y_center[tempj], kerroin, length_, temp_koko_orth);
								if (type == 2u) {
									orth_distance_precompute_3D(tempi, Nx, Nz, y_diff, x_diff, z_diff, x_center, y_center[tempj], z_center, kerroinz,
										detectors, temp_koko_orth_3D, tempk, lo, n_tempk, dec, decx);
								}
							}
							break;
						}

					}
					// The number of voxels the LOR/ray traverses
					if (type > 0u) {
						lor_orth[lo] = temp_koko_orth;
						if (type == 2u)
							lor_orth[lo + loop_var_par] = temp_koko_orth_3D;
					}
					lor[lo] = temp_koko;
				}
			}
		}
		else {

			if (fabs(y_diff) < 1e-8) {
				if (detectors.yd <= maxyy && detectors.yd >= by) {


					int32_t tempi = 0, tempk = 0, tempj = 0, iu = 0, ku = 1;
					double txu = 0., tzu = 0., tc = 0., tx0 = 0., tz0 = 0.;

					const bool skip = siddon_pre_loop_2D(bx, bz, x_diff, z_diff, maxxx, bzb, dx, dz, Nx, Nz, tempi, tempk, txu, tzu, Np, TYPE,
						detectors.zs, detectors.xs, detectors.zd, detectors.xd, tc, iu, ku, tx0, tz0);

					if (!skip) {

						const double temp_x = detectors.xs;
						detectors.xs = detectors.ys;
						detectors.ys = temp_x;

						double apu1;
						for (uint32_t ii = 0u; ii < Ny; ii++) {
							apu1 = (yy_vec[ii + 1u] - detectors.yd);
							if (apu1 > 0.) {
								tempj = static_cast<int32_t>(ii);
								break;
							}
						}
						int32_t n_tempk = tempk;
						bool z_bool = true;

						for (uint32_t ii = 0u; ii < Np; ii++) {

							temp_koko++;

							if (tx0 < tz0) {

								if (type > 0) {
									orth_distance_precompute(tempj, Ny, -x_diff, -y_diff, y_center, x_center[tempi], kerroin, length_, temp_koko_orth);
									if (type == 2u) {
										orth_distance_precompute_3D(tempj, Ny, Nz, x_diff, y_diff, z_diff, y_center, x_center[tempi], z_center, kerroinz,
											detectors, temp_koko_orth_3D, tempk, lo, n_tempk, dec, decx);
									}
								}
								tempi += iu;
								tx0 += txu;
								xyz = 1u;
								z_bool = true;
							}
							else {

								if (type > 0u) {
									orth_distance_precompute(tempj, Ny, -x_diff, -y_diff, y_center, x_center[tempi], kerroin, length_, temp_koko_orth);
									if (ii == (Np - 1u) && type == 2u) {
										orth_distance_precompute_3D(tempj, Ny, Nz, x_diff, y_diff, z_diff, y_center, x_center[tempi], z_center, kerroinz,
											detectors, temp_koko_orth_3D, tempk, lo, n_tempk, dec, decx);
									}
								}
								tempk += ku;
								tz0 += tzu;
								xyz = 3u;
								if (z_bool && tempk >= 0 && tempk < Nz) {
									n_tempk = tempk;
									z_bool = false;
								}
							}
							if (tempk < 0 || tempi < 0 || tempi >= Nx || tempk >= Nz) {
								if (type == 2u && xyz == 3u && ii != (Np - 1u)) {
									tempk -= ku;
									orth_distance_precompute_3D(tempj, Ny, Nz, x_diff, y_diff, z_diff, y_center, x_center[tempi], z_center, kerroinz,
										detectors, temp_koko_orth_3D, tempk, lo, n_tempk, dec, decx);
								}
								break;
							}

						}

						// The number of voxels the LOR/ray traverses
						if (type > 0u) {
							lor_orth[lo] = temp_koko_orth;
							if (type == 2u)
								lor_orth[lo + loop_var_par] = temp_koko_orth_3D;
						}
						lor[lo] = temp_koko;
					}

				}
				// LOR/ray doesn't traverse through the pixel space
			}
			else if (fabs(x_diff) < 1e-8) {
				if (detectors.xd <= maxxx && detectors.xd >= bx) {

					int32_t tempi = 0, tempk = 0, tempj = 0, ju = 0, ku = 1;
					double tyu = 0., tzu = 0., tc = 0., ty0 = 0., tz0 = 0.;
					const bool skip = siddon_pre_loop_2D(by, bz, y_diff, z_diff, maxyy, bzb, dy, dz, Ny, Nz, tempj, tempk, tyu, tzu, Np, TYPE,
						detectors.zs, detectors.ys, detectors.zd, detectors.yd, tc, ju, ku, ty0, tz0);

					if (!skip) {

						double apu1;
						for (uint32_t ii = 0u; ii < Nx; ii++) {
							apu1 = (xx_vec[ii + 1u] - detectors.xd);
							if (apu1 > 0.) {
								tempi = ii;
								break;
							}
						}
						int32_t n_tempk = tempk;
						bool z_bool = true;


						for (uint32_t ii = 0u; ii < Np; ii++) {

							temp_koko++;

							if (ty0 < tz0) {

								if (type > 0u) {
									orth_distance_precompute(tempi, Nx, y_diff, x_diff, x_center, y_center[tempj], kerroin, length_, temp_koko_orth);
									if (type == 2u)
										orth_distance_precompute_3D(tempi, Nx, Nz, y_diff, x_diff, z_diff, x_center, y_center[tempj], z_center, kerroinz,
											detectors, temp_koko_orth_3D, tempk, lo, n_tempk, dec, decx);
								}
								tempj += ju;
								ty0 += tyu;
								xyz = 2u;
								z_bool = true;
							}
							else {

								if (type > 0u) {
									orth_distance_precompute(tempi, Nx, y_diff, x_diff, x_center, y_center[tempj], kerroin, length_, temp_koko_orth);
									if (ii == (Np - 1u) && type == 2u) {
										orth_distance_precompute_3D(tempi, Nx, Nz, y_diff, x_diff, z_diff, x_center, y_center[tempj], z_center, kerroinz,
											detectors, temp_koko_orth_3D, tempk, lo, n_tempk, dec, decx);
									}
								}
								tempk += ku;
								tz0 += tzu;
								xyz = 3u;
								if (z_bool && tempk >= 0 && tempk < Nz) {
									n_tempk = tempk;
									z_bool = false;
								}
							}
							if (tempj < 0 || tempk < 0 || tempk >= Nz || tempj >= Ny) {
								if (type == 2u && xyz == 3u && ii != (Np - 1u)) {
									tempk -= ku;
									orth_distance_precompute_3D(tempi, Nx, Nz, y_diff, x_diff, z_diff, x_center, y_center[tempj], z_center, kerroinz,
										detectors, temp_koko_orth_3D, tempk, lo, n_tempk, dec, decx);
								}
								break;
							}

						}
						// The number of voxels the LOR/ray traverses
						if (type > 0u) {
							lor_orth[lo] = temp_koko_orth;
							if (type == 2u)
								lor_orth[lo + loop_var_par] = temp_koko_orth_3D;
						}
						lor[lo] = temp_koko;
					}
				}
				// LOR/ray doesn't traverse through the pixel space
			}
			else {

				int32_t tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 1;
				double txu = 0., tyu = 0., tzu = 0., tc = 0., tx0 = 0., ty0 = 0., tz0 = 0.;
				const bool skip = siddon_pre_loop_3D(bx, by, bz, x_diff, y_diff, z_diff, maxxx, maxyy, bzb, dx, dy, dz, Nx, Ny, Nz, tempi, tempj, tempk, tyu, txu, tzu,
					Np, TYPE, detectors, tc, iu, ju, ku, tx0, ty0, tz0);

				if (!skip) {

					int32_t n_tempk = tempk;
					bool z_bool = true;

					for (uint32_t ii = 0u; ii < Np; ii++) {

						temp_koko++;

						if (tz0 < ty0 && tz0 < tx0) {

							if (type > 0u) {
								orth_distance_precompute(tempi, Nx, y_diff, x_diff, x_center, y_center[tempj], kerroin, length_, temp_koko_orth);
								if (ii == (Np - 1u) && type == 2u) {
									orth_distance_precompute_3D(tempi, Nx, Nz, y_diff, x_diff, z_diff, x_center, y_center[tempj], z_center, kerroinz,
										detectors, temp_koko_orth_3D, tempk, lo, n_tempk, dec, decx);
								}
							}
							tempk += ku;
							tz0 += tzu;
							xyz = 3u;
							if (z_bool && tempk >= 0 && tempk < Nz) {
								n_tempk = tempk;
								z_bool = false;
							}
						}
						else if (ty0 < tx0) {

							if (type > 0u) {
								orth_distance_precompute(tempi, Nx, y_diff, x_diff, x_center, y_center[tempj], kerroin, length_, temp_koko_orth);
								if (type == 2u)
									orth_distance_precompute_3D(tempi, Nx, Nz, y_diff, x_diff, z_diff, x_center, y_center[tempj], z_center, kerroinz,
										detectors, temp_koko_orth_3D, tempk, lo, n_tempk, dec, decx);
							}
							tempj += ju;
							ty0 += tyu;
							xyz = 2u;
							z_bool = true;
						}
						else {

							if (ii == (Np - 1u) && type > 0u) {
								orth_distance_precompute(tempi, Nx, y_diff, x_diff, x_center, y_center[tempj], kerroin, length_, temp_koko_orth);
								if (type == 2u)
									orth_distance_precompute_3D(tempi, Nx, Nz, y_diff, x_diff, z_diff, x_center, y_center[tempj], z_center, kerroinz,
										detectors, temp_koko_orth_3D, tempk, lo, n_tempk, dec, decx);
							}
							else {
								tempi += iu;
								tx0 += txu;
							}
							xyz = 1u;
						}
						if (tempj < 0 || tempi < 0 || tempk < 0 || tempi >= Nx || tempj >= Ny || tempk >= Nz) {
							if (type > 0u && ii != (Np - 1u)) {
								if (xyz == 1u)
									tempi -= iu;
								else
									tempk -= ku;
								if (xyz == 1u) {
									orth_distance_precompute(tempi, Nx, y_diff, x_diff, x_center, y_center[tempj], kerroin, length_, temp_koko_orth);
									if (type == 2u) {
										orth_distance_precompute_3D(tempi, Nx, Nz, y_diff, x_diff, z_diff, x_center, y_center[tempj], z_center, kerroinz,
											detectors, temp_koko_orth_3D, tempk, lo, n_tempk, dec, decx);
									}
								}
								else if (type == 2u && xyz == 3u) {
									orth_distance_precompute_3D(tempi, Nx, Nz, y_diff, x_diff, z_diff, x_center, y_center[tempj], z_center, kerroinz,
										detectors, temp_koko_orth_3D, tempk, lo, n_tempk, dec, decx);
								}
							}
							break;
						}

					}
					// LOR/ray traverses through the pixel space
					// The number of voxels the LOR/ray traverses
					if (type > 0u) {
						lor_orth[lo] = temp_koko_orth;
						if (type == 2u)
							lor_orth[lo + loop_var_par] = temp_koko_orth_3D;
					}
					lor[lo] = temp_koko;
				}
				// LOR/ray doesn't traverse through the pixel space
			}
		}
	});
	return;
}

void mexFunction(int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray*prhs[])

{
	if (nrhs != 31)
		mexErrMsgTxt("Too few input arguments. There must be exactly 31.");

	if (nlhs != 2)
		mexErrMsgTxt("Invalid number of output arguments.  There must be exactly two.");

	const uint32_t TotSinos = (uint32_t)mxGetScalar(prhs[0]);

	const uint32_t Ny = (uint32_t)mxGetScalar(prhs[1]);

	const uint32_t Nx = (uint32_t)mxGetScalar(prhs[2]);

	const uint32_t Nz = (uint32_t)mxGetScalar(prhs[3]);

	const double d = (double)mxGetScalar(prhs[4]);

	const double dz = (double)mxGetScalar(prhs[5]);

	const double by = (double)mxGetScalar(prhs[6]);

	const double bx = (double)mxGetScalar(prhs[7]);

	const double bz = (double)mxGetScalar(prhs[8]);

	const double *z_det = (double*)mxGetData(prhs[9]);

	const double *x = (double*)mxGetData(prhs[10]);

	const double *y = (double*)mxGetData(prhs[11]);

	const double dy = (double)mxGetScalar(prhs[12]);

	const double *yy = (double*)mxGetData(prhs[13]);

	const double *xx = (double*)mxGetData(prhs[14]);

	const uint32_t NSinos = (uint32_t)mxGetScalar(prhs[15]);

	const uint32_t NSlices = (uint32_t)mxGetScalar(prhs[16]);

	const vector<double> z_det_vec(z_det, z_det + mxGetNumberOfElements(prhs[9]));

	const vector<double> yy_vec(yy, yy + mxGetNumberOfElements(prhs[13]));

	const vector<double> xx_vec(xx, xx + mxGetNumberOfElements(prhs[14]));

	const uint32_t size_x = (uint32_t)mxGetScalar(prhs[17]);

	const double zmax = (double)mxGetScalar(prhs[18]);

	const uint32_t block1 = (uint32_t)mxGetScalar(prhs[19]);

	const uint32_t blocks = (uint32_t)mxGetScalar(prhs[20]);

	const uint32_t det_per_ring = (uint32_t)mxGetScalar(prhs[21]);

	const uint16_t *L = (uint16_t*)mxGetData(prhs[22]);
	const size_t numRows = (uint32_t)mxGetM(prhs[22]);

	const uint32_t *pseudos = (uint32_t*)mxGetData(prhs[23]);
	const uint32_t pRows = (uint32_t)mxGetM(prhs[23]);

	const bool raw = (bool)mxGetScalar(prhs[24]);

	const uint32_t type = (uint32_t)mxGetScalar(prhs[25]);

	const double crystal_size = (double)mxGetScalar(prhs[26]);

	const double *x_center = (double*)mxGetData(prhs[27]);

	const double *y_center = (double*)mxGetData(prhs[28]);

	const double* z_center = (double*)mxGetData(prhs[29]);

	const double crystal_size_z = (double)mxGetScalar(prhs[30]);

	//const uint32_t projector_type = (uint32_t)mxGetScalar(prhs[25]);

	//double* x_center, *y_center;
	//double crystal_size;

	//if (projector_type == 2) {
	//	x_center = (double*)mxGetData(prhs[26]);
	//	y_center = (double*)mxGetData(prhs[27]);
	//	crystal_size = (double)mxGetScalar(prhs[28]);
	//}

	size_t loop_var_par = 1ULL;

	if (raw)
		loop_var_par = numRows / 2ULL;
	else
		loop_var_par = NSinos * size_x;

	//plhs[0] = mxCreateLogicalMatrix(loop_var_par, 1);

	//bool* discard = (bool*)mxGetLogicals(plhs[0]);

	uint16_t* lor, *lor_orth;

	if (type == 0) {

		plhs[0] = mxCreateNumericMatrix(loop_var_par, 1, mxUINT16_CLASS, mxREAL);

		lor = (uint16_t*)mxGetData(plhs[0]);

		plhs[1] = mxCreateNumericMatrix(1, 1, mxUINT16_CLASS, mxREAL);

		lor_orth = (uint16_t*)mxGetData(plhs[1]);
	}
	else {

		plhs[0] = mxCreateNumericMatrix(loop_var_par, 1, mxUINT16_CLASS, mxREAL);

		lor = (uint16_t*)mxGetData(plhs[0]);

		if (type == 2u)
			plhs[1] = mxCreateNumericMatrix(loop_var_par * 2ULL, 1, mxUINT16_CLASS, mxREAL);
		else
			plhs[1] = mxCreateNumericMatrix(loop_var_par, 1, mxUINT16_CLASS, mxREAL);

		lor_orth = (uint16_t*)mxGetData(plhs[1]);
	}

	// The maximum elements of the pixel space in both x- and y-directions
	const double maxyy = yy_vec.back();
	const double maxxx = xx_vec.back();

	improved_siddon(loop_var_par, size_x, zmax, TotSinos, lor, maxyy, maxxx, xx_vec, z_det_vec, dy, yy_vec, x, y, z_det, NSlices, Nx, Ny, Nz, 
		d, dz, bx, by, bz, block1, blocks, L, pseudos, raw, pRows, det_per_ring, type, lor_orth, crystal_size, crystal_size_z, x_center, y_center, z_center);

}