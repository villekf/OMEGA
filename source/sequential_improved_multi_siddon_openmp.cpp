#include "projector_functions.h"
#ifdef _OPENMP
#include <omp.h>
#endif

// if 0, then determines whether the LOR intercepts the FOV
constexpr int TYPE = 0;

constexpr auto THR = 0.01;

using namespace std;

void sequential_improved_siddon4(const size_t loop_var_par, const uint32_t size_x, const double zmax, double* Summ, double* rhs, const double maxyy, const double maxxx,
	const vector<double>& xx_vec, const double dy, const vector<double>& yy_vec, const double* atten, const double* x, const double* y, const double* z_det,
	const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const double dx, const double dz, const double bx, const double by, const double bz,
	const bool attenuation_correction, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t TotSinos, 
	const double epps, const double* Sino, double* osem_apu, const uint16_t *L, const uint32_t *pseudos, const size_t pRows, const uint32_t det_per_ring,
	const bool raw, const double cr_pz) {


	const uint32_t Nyx = Ny * Nx;

	const double bzb = bz + static_cast<double>(Nz) * dz;

	#pragma omp parallel for
	for (int32_t lo = 0; lo < loop_var_par; lo++) {



		Det detectors;

		int32_t tempi_a[4], tempj_a[4], tempk_a[4];
		double tx0_a[4], ty0_a[4], tz0_a[4], tc_a[4], y_diff[4], x_diff[4], z_diff[4];

		double temp = 0.;
		double ax = 0., jelppi = 0.;

		// Load the number of voxels the LOR traverses (precomputed)
		//const uint32_t Np = static_cast<uint32_t>(lor1[lo]);

		for (uint16_t lor = 0u; lor < 4u; lor++) {

			// Raw list-mode data
			if (raw) {
				//get_detector_coordinates_raw_4(det_per_ring, x, y, z_det, detectors, L, lo, pseudos, pRows);
			}
			// Sinogram data
			else {
				get_detector_coordinates_mr(x, y, z_det, size_x, detectors, xy_index, z_index, TotSinos, lo, lor + 1u, cr_pz);
			}

			// Calculate the x, y and z distances of the detector pair
			y_diff[lor] = (detectors.yd - detectors.ys);
			x_diff[lor] = (detectors.xd - detectors.xs);
			z_diff[lor] = (detectors.zd - detectors.zs);

			uint32_t Np = 0u;


			if (fabs(z_diff[lor]) < 1e-8) {

				const uint32_t tempk = z_ring(zmax, detectors.zs, static_cast<double>(NSlices));

				if (fabs(y_diff[lor]) < 1e-8) {

					if (detectors.yd <= maxyy && detectors.yd >= by) {
						uint32_t temp_ijk = 0;
						uint32_t apu = 0u;
						// Find the closest y-index value by finding the smallest y-distance between detector 2 and all the y-pixel coordinates
						for (size_t ii = 0ULL; ii < static_cast<size_t>(Ny); ii++) {
							double temp = (yy_vec[ii + 1ULL] - detectors.yd);
							if (temp > 0.) {
								apu = static_cast<uint32_t>(ii);
								break;
							}
						}
						//const double templ_ijk = d;
						temp_ijk = apu * Ny + tempk * Nx * Ny;
						temp += dx * static_cast<double>(Ny);

						if (Sino[lo] != 0.) {
							for (uint32_t k = 0; k < Nx; k++) {
								ax += (dx * osem_apu[temp_ijk + k]);
							}
						}
					}
				}
				else if (fabs(x_diff[lor]) < 1e-8) {

					if (detectors.xd <= maxxx && detectors.xd >= bx) {
						uint32_t temp_ijk = 0;

						const double element = perpendicular_elements(1, detectors.xd, xx_vec, dy, tempk, Ny, Nx, atten, attenuation_correction, temp_ijk, Nx);

						if (Sino[lo] != 0.) {
							for (uint32_t k = 0; k < Ny; k++) {
								ax += (element * osem_apu[temp_ijk + k]);
							}
							if (ax == 0.)
								ax = epps;
							const double yax = Sino[lo] / ax;
							for (uint32_t k = 0; k < Ny; k++) {
#pragma omp atomic
								rhs[tempk + k] += (element * yax);
#pragma omp atomic
								Summ[tempk + k] += element;
							}
						}
						else {
							for (uint32_t k = 0; k < Ny; k++) {
#pragma omp atomic
								Summ[tempk + k] += element;
							}
						}
					}
				}
				else {
					int32_t tempi = 0, tempj = 0, iu = 0, ju = 0;
					double txu = 0., tyu = 0., tc = 0., tx0 = 0., ty0 = 0.;

					const bool skip = siddon_pre_loop_2D(bx, by, x_diff[lor], y_diff[lor], maxxx, maxyy, dx, dy, Nx, Ny, tempi, tempj, txu, tyu, Np, TYPE,
						detectors.ys, detectors.xs, detectors.yd, detectors.xd, tc, iu, ju, tx0, ty0);

					if (tempi < 0 || tempj < 0 || tempi >= static_cast<int32_t>(Nx) || tempj >= static_cast<int32_t>(Ny))
						continue;

					const double L = sqrt(x_diff[lor] * x_diff[lor] + y_diff[lor] * y_diff[lor]);

					tempi_a[lor] = tempi, tempj_a[lor] = tempj, tempk_a[lor] = tempk;
					tx0_a[lor] = tx0, ty0_a[lor] = ty0, tc_a[lor] = tc, tz0_a[lor] = -1.;
					uint32_t tempijk = tempk * Nyx + static_cast<uint32_t>(tempj) * Nx + static_cast<uint32_t>(tempi);

					for (uint32_t ii = 0; ii < Np; ii++) {

						if (tx0 < ty0) {
							const double element = (tx0 - tc) * L;

							if (iu > 0)
								tempi++;
							else
								tempi--;
							tc = tx0;
							tx0 += txu;

							temp += element;
							ax += (element * osem_apu[temp_ijk]);
							if (attenuation_correction)
								jelppi += (element * -atten[temp_ijk]);


						}
						else {

							const double element = (ty0 - tc) * L;

							if (ju > 0)
								tempj += Nx;
							else
								tempj -= Nx;
							tc = ty0;
							ty0 += tyu;

							temp += element;
							ax += (element * osem_apu[temp_ijk]);
							if (attenuation_correction)
								jelppi += (element * -atten[temp_ijk]);
						}
					}

				}
			}
		}
		temp = 1. / temp;
		if (attenuation_correction)
			temp *= exp(jelppi);

		for (int32_t lor = 0; lor < 4; lor++) {
			//tx0 = tx0_a;
			//ty0 = ty0_a;
			//tempi = tempi_a;
			//tempj = tempj_a;
			//tc = tc_a;

			if (Sino[lo] != 0.) {
				ax *= temp;
				if (ax == 0.)
					ax = epps;
				const double yax = Sino[lo] / ax;
				for (uint32_t ii = 0; ii < Np; ii++) {

					const uint32_t temp_ijk = tempj + tempi + tempk;

					if (tx0 < ty0) {
						const double element = (tx0 - tc) * L * temp;

						if (iu > 0)
							tempi++;
						else
							tempi--;
						tc = tx0;
						tx0 += txu;


#pragma omp atomic
						rhs[temp_ijk] += (element * yax);
#pragma omp atomic
						Summ[temp_ijk] += element;


					}
					else {

						const double element = (ty0 - tc) * L * temp;

						if (ju > 0)
							tempj += Nx;
						else
							tempj -= Nx;
						tc = ty0;
						ty0 += tyu;

#pragma omp atomic
						rhs[temp_ijk] += (element * yax);
#pragma omp atomic
						Summ[temp_ijk] += element;
					}
				}
			}
			else {
				for (uint32_t ii = 0; ii < Np; ii++) {

					const uint32_t temp_ijk = tempj + tempi + tempk;

					if (tx0 < ty0) {
						const double element = (tx0 - tc) * L * temp;

						if (iu > 0)
							tempi++;
						else
							tempi--;
						tc = tx0;
						tx0 += txu;

#pragma omp atomic
						Summ[temp_ijk] += element;


					}
					else {

						const double element = (ty0 - tc) * L * temp;

						if (ju > 0)
							tempj += Nx;
						else
							tempj -= Nx;
						tc = ty0;
						ty0 += tyu;

#pragma omp atomic
						Summ[temp_ijk] += element;
					}
				}
			}
		}


					//if (attenuation_correction)
					//	att_corr_vec(elements, indices, atten, temp, static_cast<size_t>(Np));

					//for (uint32_t ii = 0; ii < Np; ii++)
					//	elements[ii] = elements[ii] * (temp);

					//pass = true;
				}
			}
			else {

				if (fabs(y_diff) < 1e-8) {
					if (detectors.yd <= maxyy && detectors.yd >= by) {

						const double apu_tz = bz - detectors.zs;
						const double apu_tx = bx - detectors.xs;
						double tx0 = (apu_tx) / (x_diff);
						double tz0 = (apu_tz) / (z_diff);
						const double txback = (maxxx - detectors.xs) / (x_diff);
						const double tzback = (bzb - detectors.zs) / (z_diff);

						const double txmin = min(tx0, txback);
						const double txmax = max(tx0, txback);
						const double tzmin = min(tz0, tzback);
						const double tzmax = max(tz0, tzback);

						const double tmin = max(txmin, tzmin);
						const double tmax = min(txmax, tzmax);

						uint32_t imin, imax, kmin, kmax;
						int iu, ku;

						if (detectors.xs < detectors.xd)
							d_g_s(tmin, txmin, tmax, txmax, imin, imax, tx0, iu, x_diff, bx, dx, detectors.xs, Nx);
						else
							s_g_d(tmin, txmin, tmax, txmax, imin, imax, tx0, iu, x_diff, bx, dx, detectors.xs, Nx);

						if (detectors.zs < detectors.zd)
							d_g_s(tmin, tzmin, tmax, tzmax, kmin, kmax, tz0, ku, z_diff, bz, dz, detectors.zs, Nz);
						else
							s_g_d(tmin, tzmin, tmax, tzmax, kmin, kmax, tz0, ku, z_diff, bz, dz, detectors.zs, Nz);

						const double L = sqrt((x_diff*x_diff + z_diff * z_diff));

						uint32_t tempi, tempj, tempk;
						double apu1;

						const double pt = ((min(tx0, tz0) + tmin) / 2.);

						tempk = Nyx * static_cast<uint32_t>(voxel_index(pt, z_diff, dz, apu_tz));
						tempi = static_cast<uint32_t>(voxel_index(pt, x_diff, dx, apu_tx));

						const double txu = dx / fabs(x_diff);
						const double tzu = dz / fabs(z_diff);

						for (uint32_t ii = 0; ii < Ny; ii++) {
							apu1 = (yy_vec[ii + 1] - detectors.yd);
							if (apu1 > 0.) {
								tempj = ii;
								break;
							}
						}

						tempj = tempj * Nx;

						double tc = tmin;

						double temp = 0.;
						uint32_t tempi_a = tempi, tempj_a = tempj, tempk_a = tempk;
						double tx0_a = tx0, tz0_a = tz0, tc_a = tc;

						for (uint32_t ii = 0; ii < Np; ii++) {

							const uint32_t temp_ijk = tempj + tempi + tempk;
							if (tx0 < tz0) {

								const double element = (tx0 - tc) * L;


								if (iu > 0)
									tempi++;
								else
									tempi--;
								tc = tx0;
								tx0 += txu;

								temp += element;
								ax += (element * osem_apu[temp_ijk]);
								if (attenuation_correction)
									jelppi += (element * -atten[temp_ijk]);
							}
							else {

								const double element = (tz0 - tc) * L;

								if (ku > 0)
									tempk += Nyx;
								else
									tempk -= Nyx;
								tc = tz0;
								tz0 += tzu;

								temp += element;
								ax += (element * osem_apu[temp_ijk]);
								if (attenuation_correction)
									jelppi += (element * -atten[temp_ijk]);
							}

						}

						temp = 1. / temp;
						tx0 = tx0_a;
						tz0 = tz0_a;
						tempi = tempi_a;
						tempk = tempk_a;
						tc = tc_a;
						if (attenuation_correction)
							temp *= exp(jelppi);


						if (Sino[lo] != 0.) {
							ax *= temp;
							if (ax == 0.)
								ax = epps;
							const double yax = Sino[lo] / ax;
							for (uint32_t ii = 0; ii < Np; ii++) {

								const uint32_t temp_ijk = tempj + tempi + tempk;
								if (tx0 < tz0) {

									const double element = (tx0 - tc) * L * temp;


									if (iu > 0)
										tempi++;
									else
										tempi--;
									tc = tx0;
									tx0 += txu;

#pragma omp atomic
									rhs[temp_ijk] += (element * yax);
#pragma omp atomic
									Summ[temp_ijk] += element;

								}
								else {

									const double element = (tz0 - tc) * L * temp;

									if (ku > 0)
										tempk += Nyx;
									else
										tempk -= Nyx;
									tc = tz0;
									tz0 += tzu;

#pragma omp atomic
									rhs[temp_ijk] += (element * yax);
#pragma omp atomic
									Summ[temp_ijk] += element;

								}

							}
						}
						else {
							for (uint32_t ii = 0; ii < Np; ii++) {

								const uint32_t temp_ijk = tempj + tempi + tempk;
								if (tx0 < tz0) {

									const double element = (tx0 - tc) * L * temp;


									if (iu > 0)
										tempi++;
									else
										tempi--;
									tc = tx0;
									tx0 += txu;

#pragma omp atomic
									Summ[temp_ijk] += element;

								}
								else {

									const double element = (tz0 - tc) * L * temp;

									if (ku > 0)
										tempk += Nyx;
									else
										tempk -= Nyx;
									tc = tz0;
									tz0 += tzu;

#pragma omp atomic
									Summ[temp_ijk] += element;

								}

							}

						}


						//if (attenuation_correction)
						//	att_corr_vec(elements, indices, atten, temp, static_cast<size_t>(Np));

						//for (uint32_t ii = 0; ii < Np; ii++)
						//	elements[ii] = (elements[ii] * temp);
						//pass = true;

					}
				}
				else if (fabs(x_diff) < 1e-8) {
					if (detectors.xd <= maxxx && detectors.xd >= bx) {

						const double apu_tz = bz - detectors.zs;
						const double apu_ty = by - detectors.ys;
						double ty0 = (apu_ty) / (y_diff);
						const double tyback = (maxyy - detectors.ys) / (y_diff);
						double tz0 = (apu_tz) / (z_diff);
						const double tzback = (bzb - detectors.zs) / (z_diff);

						const double tzmin = min(tz0, tzback);
						const double tzmax = max(tz0, tzback);
						const double tymin = min(ty0, tyback);
						const double tymax = max(ty0, tyback);

						const double tmin = max(tymin, tzmin);
						const double tmax = min(tymax, tzmax);

						uint32_t jmin, jmax, kmin, kmax;
						int ku, ju;

						if (detectors.ys < detectors.yd)
							d_g_s(tmin, tymin, tmax, tymax, jmin, jmax, ty0, ju, y_diff, by, dy, detectors.ys, Ny);
						else
							s_g_d(tmin, tymin, tmax, tymax, jmin, jmax, ty0, ju, y_diff, by, dy, detectors.ys, Ny);

						if (detectors.zs < detectors.zd)
							d_g_s(tmin, tzmin, tmax, tzmax, kmin, kmax, tz0, ku, z_diff, bz, dz, detectors.zs, Nz);
						else
							s_g_d(tmin, tzmin, tmax, tzmax, kmin, kmax, tz0, ku, z_diff, bz, dz, detectors.zs, Nz);

						const double L = sqrt((y_diff*y_diff + z_diff * z_diff));

						uint32_t tempi, tempj, tempk;
						double apu1;

						const double pt = ((min(tz0, ty0) + tmin) / 2.);

						tempk = Nyx * static_cast<uint32_t>(voxel_index(pt, z_diff, dz, apu_tz));
						tempj = Nx * static_cast<uint32_t>(voxel_index(pt, y_diff, dy, apu_ty));

						const double tzu = dz / fabs(z_diff);
						const double tyu = dy / fabs(y_diff);

						double tc = tmin;

						double temp = 0.;

						for (uint32_t ii = 0; ii < Nx; ii++) {
							apu1 = (xx_vec[ii + 1] - detectors.xd);
							if (apu1 > 0.) {
								tempi = ii;
								break;
							}
						}

						uint32_t tempi_a = tempi, tempj_a = tempj, tempk_a = tempk;
						double ty0_a = ty0, tz0_a = tz0, tc_a = tc;

						for (uint32_t ii = 0; ii < Np; ii++) {

							const uint32_t temp_ijk = tempj + tempi + tempk;
							if (ty0 < tz0) {

								const double element = (ty0 - tc) * L;


								if (ju > 0)
									tempj += Nx;
								else
									tempj -= Nx;
								tc = ty0;
								ty0 += tyu;

								temp += element;
								ax += (element * osem_apu[temp_ijk]);
								if (attenuation_correction)
									jelppi += (element * -atten[temp_ijk]);
							}
							else {

								const double element = (tz0 - tc) * L;

								if (ku > 0)
									tempk += Nyx;
								else
									tempk -= Nyx;
								tc = tz0;
								tz0 += tzu;

								temp += element;
								ax += (element * osem_apu[temp_ijk]);
								if (attenuation_correction)
									jelppi += (element * -atten[temp_ijk]);
							}

						}

						temp = 1. / temp;
						ty0 = ty0_a;
						tz0 = tz0_a;
						tempj = tempj_a;
						tempk = tempk_a;
						tc = tc_a;
						if (attenuation_correction)
							temp *= exp(jelppi);


						if (Sino[lo] != 0.) {
							ax *= temp;
							if (ax == 0.)
								ax = epps;
							const double yax = Sino[lo] / ax;
							for (uint32_t ii = 0; ii < Np; ii++) {

								const uint32_t temp_ijk = tempj + tempi + tempk;
								if (ty0 < tz0) {

									const double element = (ty0 - tc) * L * temp;


									if (ju > 0)
										tempj += Nx;
									else
										tempj -= Nx;
									tc = ty0;
									ty0 += tyu;

#pragma omp atomic
									rhs[temp_ijk] += (element * yax);
#pragma omp atomic
									Summ[temp_ijk] += element;

								}
								else {

									const double element = (tz0 - tc) * L * temp;

									if (ku > 0)
										tempk += Nyx;
									else
										tempk -= Nyx;
									tc = tz0;
									tz0 += tzu;

#pragma omp atomic
									rhs[temp_ijk] += (element * yax);
#pragma omp atomic
									Summ[temp_ijk] += element;

								}

							}
						}
						else {
							for (uint32_t ii = 0; ii < Np; ii++) {

								const uint32_t temp_ijk = tempj + tempi + tempk;
								if (ty0 < tz0) {

									const double element = (ty0 - tc) * L * temp;


									if (ju > 0)
										tempj += Nx;
									else
										tempj -= Nx;
									tc = ty0;
									ty0 += tyu;
#pragma omp atomic
									Summ[temp_ijk] += element;

								}
								else {

									const double element = (tz0 - tc) * L * temp;

									if (ku > 0)
										tempk += Nyx;
									else
										tempk -= Nyx;
									tc = tz0;
									tz0 += tzu;

#pragma omp atomic
									Summ[temp_ijk] += element;

								}

							}

						}

						//if (attenuation_correction)
						//	att_corr_vec(elements, indices, atten, temp, static_cast<size_t>(Np));

					}
				}
				else {

					const double apu_tz = bz - detectors.zs;
					const double apu_ty = by - detectors.ys;
					const double apu_tx = bx - detectors.xs;
					double tx0 = (apu_tx) / (x_diff);
					double ty0 = (apu_ty) / (y_diff);
					double tz0 = (apu_tz) / (z_diff);
					const double txback = (maxxx - detectors.xs) / (x_diff);
					const double tyback = (maxyy - detectors.ys) / (y_diff);
					const double tzback = (bzb - detectors.zs) / (z_diff);

					const double txmin = min(tx0, txback);
					const double txmax = max(tx0, txback);
					const double tymin = min(ty0, tyback);
					const double tymax = max(ty0, tyback);
					const double tzmin = min(tz0, tzback);
					const double tzmax = max(tz0, tzback);

					const double tmin = max(max(txmin, tzmin), tymin);
					const double tmax = min(min(txmax, tzmax), tymax);

					uint32_t imin, imax, jmin, jmax, kmin, kmax;
					int iu, ju, ku;


					if (detectors.xs < detectors.xd)
						d_g_s(tmin, txmin, tmax, txmax, imin, imax, tx0, iu, x_diff, bx, dx, detectors.xs, Nx);
					else
						s_g_d(tmin, txmin, tmax, txmax, imin, imax, tx0, iu, x_diff, bx, dx, detectors.xs, Nx);

					if (detectors.ys < detectors.yd)
						d_g_s(tmin, tymin, tmax, tymax, jmin, jmax, ty0, ju, y_diff, by, dy, detectors.ys, Ny);
					else
						s_g_d(tmin, tymin, tmax, tymax, jmin, jmax, ty0, ju, y_diff, by, dy, detectors.ys, Ny);

					if (detectors.zs < detectors.zd)
						d_g_s(tmin, tzmin, tmax, tzmax, kmin, kmax, tz0, ku, z_diff, bz, dz, detectors.zs, Nz);
					else
						s_g_d(tmin, tzmin, tmax, tzmax, kmin, kmax, tz0, ku, z_diff, bz, dz, detectors.zs, Nz);

					const double L = sqrt(x_diff*x_diff + z_diff * z_diff + y_diff * y_diff);

					uint32_t tempi, tempj, tempk;

					const double pt = ((min(min(tz0, ty0), tx0) + tmin) / 2.);

					tempk = Nyx * static_cast<uint32_t>(voxel_index(pt, z_diff, dz, apu_tz));
					tempj = Nx * static_cast<uint32_t>(voxel_index(pt, y_diff, dy, apu_ty));
					tempi = static_cast<uint32_t>(voxel_index(pt, x_diff, dx, apu_tx));

					const double tzu = dz / fabs(z_diff);
					const double tyu = dy / fabs(y_diff);
					const double txu = dx / fabs(x_diff);

					double tc = tmin;

					double temp = 0.;

					uint32_t tempi_a = tempi, tempj_a = tempj, tempk_a = tempk;
					double ty0_a = ty0, tz0_a = tz0, tc_a = tc, tx0_a = tx0;

					for (uint32_t ii = 0; ii < Np; ii++) {
						const uint32_t temp_ijk = tempj + tempi + tempk;
						if (tz0 < ty0 && tz0 < tx0) {

							const double element = (tz0 - tc) * L;

							if (ku > 0)
								tempk += Nyx;
							else
								tempk -= Nyx;
							tc = tz0;
							tz0 += tzu;

							temp += element;
							ax += (element * osem_apu[temp_ijk]);
							if (attenuation_correction)
								jelppi += (element * -atten[temp_ijk]);
						}
						else if (ty0 < tx0) {
							const double element = (ty0 - tc) * L;

							if (ju > 0)
								tempj += Nx;
							else
								tempj -= Nx;
							tc = ty0;
							ty0 += tyu;

							temp += element;
							ax += (element * osem_apu[temp_ijk]);
							if (attenuation_correction)
								jelppi += (element * -atten[temp_ijk]);
						}
						else {
							const double element = (tx0 - tc) * L;

							if (iu > 0)
								tempi++;
							else
								tempi--;
							tc = tx0;
							tx0 += txu;

							temp += element;
							ax += (element * osem_apu[temp_ijk]);
							if (attenuation_correction)
								jelppi += (element * -atten[temp_ijk]);
						}

					}

					temp = 1. / temp;
					ty0 = ty0_a;
					tx0 = tx0_a;
					tz0 = tz0_a;
					tempi = tempi_a;
					tempj = tempj_a;
					tempk = tempk_a;
					tc = tc_a;
					if (attenuation_correction)
						temp *= exp(jelppi);


					if (Sino[lo] != 0.) {
						ax *= temp;
						if (ax == 0.)
							ax = epps;
						const double yax = Sino[lo] / ax;
						for (uint32_t ii = 0; ii < Np; ii++) {
							const uint32_t temp_ijk = tempj + tempi + tempk;
							if (tz0 < ty0 && tz0 < tx0) {

								const double element = (tz0 - tc) * L * temp;

								if (ku > 0)
									tempk += Nyx;
								else
									tempk -= Nyx;
								tc = tz0;
								tz0 += tzu;

#pragma omp atomic
								rhs[temp_ijk] += (element * yax);
#pragma omp atomic
								Summ[temp_ijk] += element;
							}
							else if (ty0 < tx0) {
								const double element = (ty0 - tc) * L * temp;

								if (ju > 0)
									tempj += Nx;
								else
									tempj -= Nx;
								tc = ty0;
								ty0 += tyu;

#pragma omp atomic
								rhs[temp_ijk] += (element * yax);
#pragma omp atomic
								Summ[temp_ijk] += element;
							}
							else {
								const double element = (tx0 - tc) * L * temp;

								if (iu > 0)
									tempi++;
								else
									tempi--;
								tc = tx0;
								tx0 += txu;

#pragma omp atomic
								rhs[temp_ijk] += (element * yax);
#pragma omp atomic
								Summ[temp_ijk] += element;
							}

						}

					}
					else {
						for (uint32_t ii = 0; ii < Np; ii++) {
							const uint32_t temp_ijk = tempj + tempi + tempk;
							if (tz0 < ty0 && tz0 < tx0) {

								const double element = (tz0 - tc) * L * temp;

								if (ku > 0)
									tempk += Nyx;
								else
									tempk -= Nyx;
								tc = tz0;
								tz0 += tzu;

#pragma omp atomic
								Summ[temp_ijk] += element;
							}
							else if (ty0 < tx0) {
								const double element = (ty0 - tc) * L * temp;

								if (ju > 0)
									tempj += Nx;
								else
									tempj -= Nx;
								tc = ty0;
								ty0 += tyu;

#pragma omp atomic
								Summ[temp_ijk] += element;
							}
							else {
								const double element = (tx0 - tc) * L * temp;

								if (iu > 0)
									tempi++;
								else
									tempi--;
								tc = tx0;
								tx0 += txu;

#pragma omp atomic
								Summ[temp_ijk] += element;
							}

						}

					}

					//if (attenuation_correction)
					//	att_corr_vec(elements, indices, atten, temp, static_cast<size_t>(Np));
				}
			}
		}
	}
}