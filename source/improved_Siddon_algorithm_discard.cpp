/**************************************************************************
* This function is used to check for the number of voxels each LOR/ray
* traverses and also whether the LOR/ray actually intersects with the pixel
* space.
* Raw list-mode data and sinogram data have slightly different versions.
* Output is the number of voxels the ray has traversed (if the LOR does not
* traverse the pixel space, this value will be 0).
*
* Uses OpenMP for parallelization.
*
* Copyright (C) 2020 Ville-Veikko Wettenhovi
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

const static int TYPE = 0;

// Whether to use the OpenMP code or not
const static bool OMP = false;

// Using non-OpenMP with either precomputation or without
const static bool PRECOMPUTE = false;


const static bool DISCARD = true;


using namespace std;

void improved_siddon_precomputation_phase(const int64_t loop_var_par, const uint32_t size_x, const double zmax, const uint32_t TotSinos, uint16_t* lor, const double maxyy,
	const double maxxx, const std::vector<double>& xx_vec, const std::vector<double>& z_det_vec, const double dy, const std::vector<double>& yy_vec,
	const double* x, const double* y, const double* z_det, const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const double dx, const double dz,
	const double bx, const double by, const double bz, const uint32_t block1, const uint32_t blocks, const uint16_t* L, const uint32_t* pseudos,
	const bool raw, const uint32_t pRows, const uint32_t det_per_ring, const uint32_t type, uint16_t* lor_orth, uint16_t* lor_vol, const double crystal_size, const double crystal_size_z,
	const double* x_center, const double* y_center, const double* z_center, const double bmin, const double bmax, const double Vmax, const double* V, const uint32_t nCores) {

	if (nCores == 1U)
		setThreads();
	else
		omp_set_num_threads(nCores);

	const double bzb = bz + static_cast<double>(Nz) * dz;

	const double local_sino = 0.;
	const double* osem_apu = nullptr;
	double ax = 0.;
	const bool no_norm = false;
	vector<double> elements;
	vector<uint32_t> v_indices;
	double* Summ = nullptr, *rhs = nullptr;
	const uint32_t Nyx = Nx * Ny;
	std::vector<double> store_elements;
	std::vector<uint32_t> store_indices;
	const bool RHS = false, SUMMA = false;
	const uint32_t tid = 0U;
	size_t* indi = 0;

#pragma omp parallel for ordered schedule(dynamic)
	for (int64_t lo = 0LL; lo < loop_var_par; lo++) {
		Det detectors;
		uint32_t ind = 0U;

		if (raw) {
			const uint32_t detektorit1 = static_cast<uint32_t>(L[lo * 2ULL]) - 1u;
			const uint32_t detektorit2 = static_cast<uint32_t>(L[lo * 2ULL + 1ULL]) - 1u;

			if (detektorit1 == detektorit2)
				continue;

			const uint32_t loop1 = ((detektorit1) / det_per_ring);
			const uint32_t loop2 = ((detektorit2) / det_per_ring);

			if (loop1 == loop2) {
				if (loop1 > blocks || loop1 < block1 || loop2 > blocks || loop2 < block1)
					continue;
				detectors.zs = z_det_vec[loop1];
				detectors.zd = detectors.zs;
			}
			else {
				if ((loop1 > blocks && loop2 > blocks) || (loop1 < block1 && loop2 < block1))
					continue;
				detectors.zs = z_det_vec[loop1];
				detectors.zd = z_det_vec[loop2];
			}

			detectors.xs = x[detektorit1 - det_per_ring * (loop1)];
			detectors.xd = x[detektorit2 - det_per_ring * (loop2)];
			detectors.ys = y[detektorit1 - det_per_ring * (loop1)];
			detectors.yd = y[detektorit2 - det_per_ring * (loop2)];
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


		const double x_diff = (detectors.xd - detectors.xs);
		const double y_diff = (detectors.yd - detectors.ys);
		const double z_diff = (detectors.zd - detectors.zs);
		if ((y_diff == 0. && x_diff == 0. && z_diff == 0.) || (y_diff == 0. && x_diff == 0.))
			continue;

		double kerroinz = 0., kerroin = 0.;
		double temp = 0.;

		size_t temp_koko_orth = 0ULL;
		uint16_t temp_koko = 0u;
		size_t temp_koko_vol = 0ULL;
		size_t temp_koko_orth_3D = 0ULL;

		uint32_t Np = 0u;
		uint8_t xyz = 0u;

		if (type == 1u || type == 2u) {
			if (type == 2u) {
				kerroinz = norm(x_diff, y_diff, z_diff) * crystal_size_z;
			}
			kerroin = norm(x_diff, y_diff, z_diff) * crystal_size;
		}
		else if (type == 3u)
			kerroin = norm(x_diff, y_diff, z_diff);

		if (fabs(z_diff) < 1e-8 && (fabs(y_diff) < 1e-8 || fabs(x_diff) < 1e-8)) {

			const uint32_t tempk = z_ring(zmax, detectors.zs, static_cast<double>(NSlices));
			if (tempk >= Nz)
				continue;

			if (fabs(y_diff) < 1e-8) {


				if (detectors.yd <= maxyy && detectors.yd >= by) {
					// The number of voxels the LOR/ray traverses

					if (type == 1u || type == 2u) {
						double temppi = detectors.xs;
						detectors.xs = detectors.ys;
						detectors.ys = temppi;
						orth_perpendicular_precompute(Nx, Ny, detectors.yd, yy_vec, y_center, x_center[0], z_center, kerroin, temp_koko_orth,
							detectors, y_diff, x_diff, z_diff, tempk);
						lor_orth[lo] = static_cast<uint16_t>(temp_koko_orth);
						if (type == 2u) {
							orth_perpendicular_precompute_3D(Nx, Ny, Nz, detectors.yd, yy_vec, y_center, x_center[0], z_center, kerroinz, temp_koko_orth_3D,
								detectors, y_diff, x_diff, z_diff, tempk);
							lor_orth[lo + loop_var_par] = static_cast<uint16_t>(temp_koko_orth_3D);
						}
					}
					else if (type == 3u) {
						double temppi = detectors.xs;
						detectors.xs = detectors.ys;
						detectors.ys = temppi;
						volume_perpendicular_precompute(Nx, Ny, Nz, detectors.yd, yy_vec, y_center, x_center[0], z_center, kerroin, temp_koko_vol,
							detectors, y_diff, x_diff, z_diff, tempk, bmin, bmax, Vmax, V);
						lor_vol[lo] = static_cast<uint16_t>(temp_koko_vol);
					}
					lor[lo] = static_cast<uint16_t>(Nx);
				}
				// LOR/ray doesn't traverse through the pixel space
			}
			else if (fabs(x_diff) < 1e-8) {

				if (detectors.xd <= maxxx && detectors.xd >= bx) {
					// The number of voxels the LOR/ray traverses
					if (type == 1u || type == 2u) {
						orth_perpendicular_precompute(Ny, Nx, detectors.xd, xx_vec, x_center, y_center[0], z_center, kerroin, temp_koko_orth, detectors,
							x_diff, y_diff, z_diff, tempk);
						lor_orth[lo] = static_cast<uint16_t>(temp_koko_orth);
						if (type == 2u) {
							orth_perpendicular_precompute_3D(Ny, Nx, Nz, detectors.xd, xx_vec, x_center, y_center[0], z_center, kerroinz, temp_koko_orth_3D,
								detectors, x_diff, y_diff, z_diff, tempk);
							lor_orth[lo + loop_var_par] = static_cast<uint16_t>(temp_koko_orth_3D);
						}
					}
					else if (type == 3u) {
						volume_perpendicular_precompute(Ny, Nx, Nz, detectors.xd, xx_vec, x_center, y_center[0], z_center, kerroin, temp_koko_vol,
							detectors, x_diff, y_diff, z_diff, tempk, bmin, bmax, Vmax, V);
						lor_vol[lo] = static_cast<uint16_t>(temp_koko_vol);
					}
					lor[lo] = static_cast<uint16_t>(Ny);
				}
				// LOR/ray doesn't traverse through the pixel space
			}
		}
		else {
			int32_t tempi = 0, tempj = 0, tempk = 0, iu = 0, ju = 0, ku = 0;
			double txu = 0., tyu = 0., tzu = 0., tc = 0., tx0 = 1e8, ty0 = 1e8, tz0 = 1e8;
			bool skip = false;

			if (std::fabs(z_diff) < 1e-8) {
				tempk = z_ring(zmax, detectors.zs, static_cast<double>(NSlices));
				skip = siddon_pre_loop_2D(bx, by, x_diff, y_diff, maxxx, maxyy, dx, dy, Nx, Ny, tempi, tempj, txu, tyu, Np, TYPE,
					detectors.ys, detectors.xs, detectors.yd, detectors.xd, tc, iu, ju, tx0, ty0);
			}
			//Detectors on different rings (e.g. oblique sinograms)
			else if (std::fabs(y_diff) < 1e-8) {
				skip = siddon_pre_loop_2D(bx, bz, x_diff, z_diff, maxxx, bzb, dx, dz, Nx, Nz, tempi, tempk, txu, tzu, Np, TYPE,
					detectors.zs, detectors.xs, detectors.zd, detectors.xd, tc, iu, ku, tx0, tz0);
				if (detectors.yd > maxyy || detectors.yd < by)
					skip = true;
				tempj = perpendicular_start(by, detectors.yd, dy, Ny);
			}
			else if (std::fabs(x_diff) < 1e-8) {
				skip = siddon_pre_loop_2D(by, bz, y_diff, z_diff, maxyy, bzb, dy, dz, Ny, Nz, tempj, tempk, tyu, tzu, Np, TYPE,
					detectors.zs, detectors.ys, detectors.zd, detectors.yd, tc, ju, ku, ty0, tz0);
				if (detectors.xd > maxxx || detectors.xd < bx)
					skip = true;
				tempi = perpendicular_start(bx, detectors.xd, dx, Nx);
			}
			else {
				skip = siddon_pre_loop_3D(bx, by, bz, x_diff, y_diff, z_diff, maxxx, maxyy, bzb, dx, dy, dz, Nx, Ny, Nz, tempi, tempj, tempk, tyu, txu, tzu,
					Np, TYPE, detectors, tc, iu, ju, ku, tx0, ty0, tz0);
			}

			if (!skip) {

				int alku, loppu;
				if (type > 0u) {
					if (type == 1u || type == 2u) {
						alku = tempk + 1;
						loppu = tempk;
						orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
							tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
							PRECOMPUTE, DISCARD, rhs, Summ, indi, elements, v_indices, temp_koko_orth, Ny, Nx, alku, iu, ju, loppu, store_elements, store_indices, tid, ind);
					}
					if (type > 1u) {
						alku = Nz;
						loppu = 0;
						if (ku > 0) {
							alku = tempk + 1;
						}
						else if (ku < 0) {
							loppu = tempk;
						}
						if (type == 2u) {
							orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
								tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroinz, no_norm, RHS, SUMMA, OMP,
								PRECOMPUTE, DISCARD, rhs, Summ, indi, elements, v_indices, temp_koko_orth_3D, Ny, Nx, alku, iu, ju, loppu, store_elements, store_indices, tid, ind);
						}
						if (type == 3u) {
							volume_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
								tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
								PRECOMPUTE, DISCARD, rhs, Summ, indi, elements, v_indices, temp_koko_vol, Ny, Nx, alku, iu, ju, loppu, store_elements, store_indices, tid, ind, bmax, bmin, Vmax, V);
						}
					}
				}

				for (uint32_t ii = 0u; ii < Np; ii++) {

					temp_koko++;

					if (tx0 < ty0 && tx0 < tz0) {
						tempi += iu;
						tx0 += txu;
						xyz = 1u;
					}
					else if (ty0 < tz0) {
						tempj += ju;
						ty0 += tyu;
						xyz = 2u;
					}
					else {
						tempk += ku;
						tz0 += tzu;
						xyz = 3u;
						if (type == 1u || type == 2u) {
							if (tempk < Nz && tempk >= 0) {
								alku = tempk + 1;
								loppu = tempk;
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
									tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
									PRECOMPUTE, DISCARD, rhs, Summ, indi, elements, v_indices, temp_koko_orth, Ny, Nx, alku, iu, ju, loppu, store_elements, store_indices, tid, ind);
							}
						}
						if (type > 1u) {
							if (tempk < Nz && tempk >= 0) {
								alku = tempk + 1;
								loppu = tempk;
								if (type == 2u) {
									orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
										tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroinz, no_norm, RHS, SUMMA, OMP,
										PRECOMPUTE, DISCARD, rhs, Summ, indi, elements, v_indices, temp_koko_orth_3D, Ny, Nx, alku, iu, ju, loppu, store_elements, store_indices, tid, ind);
								}
								if (type == 3u) {
									volume_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
										tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
										PRECOMPUTE, DISCARD, rhs, Summ, indi, elements, v_indices, temp_koko_vol, Ny, Nx, alku, iu, ju, loppu, store_elements, store_indices, tid, ind, bmax, bmin, Vmax, V);
								}
							}
						}
					}
					if (tempj < 0 || tempi < 0 || tempk < 0 || tempi >= static_cast<int32_t>(Nx) || tempj >= static_cast<int32_t>(Ny) || tempk >= static_cast<int32_t>(Nz)) {
						if (xyz < 3 && type > 1u) {
							if (xyz == 1)
								tempi -= iu;
							else if (xyz == 2)
								tempj -= ju;
							if ((tempk >= (Nz - 1) && ku > 0) || (tempk <= 0 && ku < 0))
								break;
							else
								tempk += ku;
							alku = Nz;
							loppu = 0;
							if (ku > 0) {
								loppu = tempk;
							}
							else if (ku < 0) {
								alku = tempk + 1;
							}
							if (type == 2u) {
								orth_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
									tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroinz, no_norm, RHS, SUMMA, OMP,
									PRECOMPUTE, DISCARD, rhs, Summ, indi, elements, v_indices, temp_koko_orth_3D, Ny, Nx, alku, iu, ju, loppu, store_elements, store_indices, tid, ind);
							}
							if (type == 3u) {
								volume_distance_3D_full(tempi, Nx, Nz, y_diff, x_diff, z_diff, y_center, x_center, z_center, temp, 1u,
									tempj, tempk, local_sino, ax, osem_apu, detectors, Nyx, kerroin, no_norm, RHS, SUMMA, OMP,
									PRECOMPUTE, DISCARD, rhs, Summ, indi, elements, v_indices, temp_koko_vol, Ny, Nx, alku, iu, ju, loppu, store_elements, store_indices, tid, ind, bmax, bmin, Vmax, V);
							}
						}
						break;
					}
				}
				// The number of voxels the LOR/ray traverses
				if (type == 1u || type == 2u) {
					lor_orth[lo] = static_cast<uint16_t>(temp_koko_orth);
					if (type == 2u)
						lor_orth[lo + loop_var_par] = static_cast<uint16_t>(temp_koko_orth_3D);
				}
				else if (type == 3u)
					lor_vol[lo] = static_cast<uint16_t>(temp_koko_vol);
				lor[lo] = temp_koko;
			}
		}
	}
	
}