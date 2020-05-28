/**************************************************************************
* Includes various functions required by the original Siddon, improved 
* Siddon and orthogonal distance based ray tracers.
* 
* References: 
* [1] Jacobs, F., Sundermann, E., De Sutter, B., Christiaens, M. Lemahieu, 
* I. (1998). A Fast Algorithm to Calculate the Exact Radiological Path 
* CCough a Pixel or Voxel Space. Journal of computing and information 
* technology, 6 (1), 89-94.
*
* Copyright (C) 2019 Ville-Veikko Wettenhovi
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
// Normalized distances below this are discarded in orthogonal ray tracer
const static auto CC = 1e3;


// Compute the orthogonal distance between a point and a line (in 3D)
double compute_element_volume_3D(Det detectors, const double xl, const double yl, const double zl, const double kerroin,
	const double xp, const double yp, const double zp) {

	double x1, y1, z1, x0, y0, z0;

	x0 = xp - detectors.xs;
	y0 = yp - detectors.ys;
	z0 = zp - detectors.zs;
	//x1 = xp - detectors.xd;
	//y1 = yp - detectors.yd;
	//z1 = zp - detectors.zd;

	// Cross product
	x1 = yl * z0 - zl * y0;
	y1 = zl * x0 - xl * z0;
	z1 = xl * y0 - yl * x0;
	//x1 = y1 * z0 - z1 * y0;
	//y1 = z1 * x0 - x1 * z0;
	//z1 = x1 * y0 - y1 * x0;

	// Distance
	return (norm(x1, y1, z1) / kerroin);
}

// compute the orthogonal distance, perpendicular detectors
void volume_perpendicular_3D(const double dd, const std::vector<double> vec, const uint32_t z_ring, const uint32_t N1, const uint32_t N2, 
	const uint32_t Nz, const uint32_t Nyx, const uint32_t d_N, const uint32_t d_NN, const Det detectors, const double xl, const double yl, 
	const double zl, const double* center1, const double center2, const double* z_center, const double crystal_size_z, int& hpk, double& temp, 
	uint32_t& tempk, mwIndex* indices, double* elements, const uint64_t Np) {
	uint32_t apu = 0u;
	// Find the closest y-index value by finding the smallest y-distance between detector 2 and all the y-pixel coordinates
	for (size_t ii = 0ULL; ii < static_cast<size_t>(N2); ii++) {
		double temp = (vec[ii + 1ULL] - dd);
		if (temp > 0.) {
			apu = static_cast<uint32_t>(ii);
			break;
		}
	}
	tempk = apu * d_N + z_ring * Nyx;
	double jelppi = 0.;
	for (int32_t zz = static_cast<int32_t>(z_ring); zz >= 0; zz--) {
		for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
			double d_ort = compute_element_volume_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (d_ort <= CC)
				break;
			uint32_t local_ind = uu * d_N + zz * Nyx;
			for (uint32_t kk = 0u; kk < N1; kk++) {
				indices[Np + hpk] = static_cast<mwIndex>(static_cast<int64_t>(local_ind) + static_cast<int64_t>(kk) * static_cast<int64_t>(d_NN));
				elements[Np + hpk] = d_ort;
				hpk++;
			}
			temp += d_ort * static_cast<double>(N2);
		}
		for (uint32_t uu = apu + 1u; uu < N2; uu++) {
			double d_ort = compute_element_volume_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (d_ort <= CC)
				break;
			uint32_t local_ind = uu * d_N + zz * Nyx;
			temp += d_ort;
			for (uint32_t kk = 0u; kk < N1; kk++) {
				indices[Np + hpk] = static_cast<mwIndex>(static_cast<int64_t>(local_ind) + static_cast<int64_t>(kk) * static_cast<int64_t>(d_NN));
				elements[Np + hpk] = d_ort;
				hpk++;
			}
			temp += d_ort * static_cast<double>(N2);
		}
	}
	for (uint32_t zz = z_ring + 1u; zz < Nz; zz++) {
		for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
			double d_ort = compute_element_volume_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (d_ort <= CC)
				break;
			uint32_t local_ind = uu * d_N + zz * Nyx;
			temp += d_ort;
			for (uint32_t kk = 0u; kk < N1; kk++) {
				indices[Np + hpk] = static_cast<mwIndex>(static_cast<int64_t>(local_ind) + static_cast<int64_t>(kk) * static_cast<int64_t>(d_NN));
				elements[Np + hpk] = d_ort;
				hpk++;
			}
			temp += d_ort * static_cast<double>(N2);
		}
		for (uint32_t uu = apu + 1; uu < N2; uu++) {
			double d_ort = compute_element_volume_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (d_ort <= CC)
				break;
			uint32_t local_ind = uu * d_N + zz * Nyx;
			temp += d_ort;
			for (uint32_t kk = 0u; kk < N1; kk++) {
				indices[Np + hpk] = static_cast<mwIndex>(static_cast<int64_t>(local_ind) + static_cast<int64_t>(kk) * static_cast<int64_t>(d_NN));
				elements[Np + hpk] = d_ort;
				hpk++;
			}
			temp += d_ort * static_cast<double>(N2);
		}
	}
	// Probability
	temp = 1. / temp;
}


// Compute the total distance (and optionally forward projection) for the orthogonal ray (2D)
void volume_distance_3D_full(int32_t tempi, const uint32_t Nx, const uint32_t Nz, const double y_diff, const double x_diff, const double z_diff,
	const double* y_center, const double* x_center, const double* z_center, double& temp, const uint32_t NN, int32_t tempj, int32_t tempk, 
	const double local_sino, double& ax, const double* osem_apu, const Det detectors, const uint32_t Nyx, const double kerroin, 
	const bool no_norm, const bool RHS, const bool SUMMA, const bool OMP, const bool PRECOMP, const bool DISCARD, double* rhs, double* Summ, mwIndex* indices,
	std::vector<double>& elements, std::vector<uint32_t>& v_indices, size_t& idx, const uint32_t Ny, const uint32_t N1, const int start,
	const int32_t iu, const int32_t ju, const int loppu, std::vector<double>& store_elements, std::vector<uint32_t>& store_indices,
	const uint32_t tid, uint32_t& ind, const double bmax, const double bmin, const double Vmax, const double* V, uint64_t N2, uint64_t N22) {

	if (RHS || SUMMA) {
		for (int32_t uu = 0; uu < ind; uu++) {
			double local_ele = store_elements[tid + uu];
			uint32_t local_ind = store_indices[tid + uu];
			computeIndices(RHS, SUMMA, OMP, PRECOMP, DISCARD, local_ele, temp, ax, no_norm, Summ, rhs,
				local_sino, osem_apu, N2, indices, elements, v_indices, idx, local_ind, N22);
		}
	}
	else {
		int yy1 = 0;
		int uu1 = 0;
		int yy2 = 0;
		int uu2 = 0;
		int alku_y1 = tempj;
		int alku_x1 = tempi;
		int alku_y2 = tempj - 1;
		int alku_x2 = tempi - 1;
		bool breikki1 = false;
		bool breikki2 = false;
		bool breikki3 = false;
		bool breikki4 = false;
		for (int zz = tempk; zz < start; zz++) {
			yy1 = 0;
			uu1 = 0;
			yy2 = 0;
			uu2 = 0;
			alku_y1 = tempj;
			alku_x1 = tempi;
			alku_y2 = tempj - 1;
			alku_x2 = tempi - 1;
			breikki1 = false;
			breikki2 = false;
			breikki3 = false;
			breikki4 = false;
			for (yy1 = alku_y1; yy1 < Ny; yy1++) {
				int xx = 0;
				int incr = 0;
				float prev_local = 1.f;
				for (xx = alku_x1; xx < Nx; xx++) {
					double local_ele = compute_element_volume_3D(detectors, x_diff, y_diff, z_diff, kerroin, x_center[xx], y_center[yy1], z_center[zz]);
					if (local_ele >= bmax && incr > 0 && prev_local < local_ele) {
						if (xx == alku_x1 + 1) {
							breikki1 = true;
						}
						break;
					}
					else if (local_ele >= bmax) {
						incr++;
						prev_local = local_ele;
						continue;
					}
					incr = 1;
					prev_local = local_ele;
					if (local_ele < bmin)
						local_ele = Vmax;
					else
						local_ele = V[static_cast<size_t>(std::round((local_ele - bmin) * CC))];
					const uint32_t local_ind = compute_ind_orth_mfree_3D(static_cast<uint32_t>(xx), yy1 * N1, static_cast<uint32_t>(zz), NN, Nyx);
					computeIndices(RHS, SUMMA, OMP, PRECOMP, DISCARD, local_ele, temp, ax, no_norm, Summ, rhs,
						local_sino, osem_apu, N2, indices, elements, v_indices, idx, local_ind, N22);
					if (!DISCARD && !PRECOMP) {
						store_elements[tid + ind] = local_ele;
						store_indices[tid + ind] = local_ind;
						ind++;
					}
				}
				uu1 = xx;
				xx = 0;
				incr = 0;
				prev_local = 1.f;
				for (xx = alku_x2; xx >= 0; xx--) {
					double local_ele = compute_element_volume_3D(detectors, x_diff, y_diff, z_diff, kerroin, x_center[xx], y_center[yy1], z_center[zz]);
					if (local_ele >= bmax && incr > 0 && prev_local < local_ele) {
						if (xx == alku_x2 - 1) {
							breikki2 = true;
						}
						break;
					}
					else if (local_ele >= bmax) {
						incr++;
						prev_local = local_ele;
						continue;
					}
					incr = 1;
					prev_local = local_ele;
					if (local_ele < bmin)
						local_ele = Vmax;
					else
						local_ele = V[static_cast<size_t>(std::round((local_ele - bmin) * CC))];
					const uint32_t local_ind = compute_ind_orth_mfree_3D(static_cast<uint32_t>(xx), yy1 * N1, static_cast<uint32_t>(zz), NN, Nyx);
					computeIndices(RHS, SUMMA, OMP, PRECOMP, DISCARD, local_ele, temp, ax, no_norm, Summ, rhs,
						local_sino, osem_apu, N2, indices, elements, v_indices, idx, local_ind, N22);
					if (!DISCARD && !PRECOMP) {
						store_elements[tid + ind] = local_ele;
						store_indices[tid + ind] = local_ind;
						ind++;
					}
				}
				uu2 = xx;
				if (iu > 0) {
					if (ju > 0) {
						alku_x1 = uu1 - 1;
						alku_x2 = uu1 - 2;
					}
					else {
						alku_x2 = uu2 + 1;
						alku_x1 = uu2 + 2;
					}
				}
				else {
					if (ju > 0) {
						alku_x2 = uu2 + 1;
						alku_x1 = uu2 + 2;
					}
					else {
						alku_x1 = uu1 - 1;
						alku_x2 = uu1 - 2;
					}
				}
				if (breikki1 && breikki2) {
					if (yy1 == alku_y1)
						breikki3 = true;
					break;
				}
			}
			breikki1 = false;
			breikki2 = false;
			alku_x1 = tempi;
			alku_x2 = tempi - 1;
			for (yy2 = alku_y2; yy2 >= 0; yy2--) {
				int xx = 0;
				int incr = 0;
				float prev_local = 1.f;
				for (xx = alku_x1; xx < Nx; xx++) {
					double local_ele = compute_element_volume_3D(detectors, x_diff, y_diff, z_diff, kerroin, x_center[xx], y_center[yy2], z_center[zz]);
					if (local_ele >= bmax && incr > 0 && prev_local < local_ele) {
						if (xx == alku_x1 + 1) {
							breikki1 = true;
						}
						break;
					}
					else if (local_ele >= bmax) {
						incr++;
						prev_local = local_ele;
						continue;
					}
					incr = 1;
					prev_local = local_ele;
					if (local_ele < bmin)
						local_ele = Vmax;
					else
						local_ele = V[static_cast<size_t>(std::round((local_ele - bmin) * CC))];
					const uint32_t local_ind = compute_ind_orth_mfree_3D(static_cast<uint32_t>(xx), yy2 * N1, static_cast<uint32_t>(zz), NN, Nyx);
					computeIndices(RHS, SUMMA, OMP, PRECOMP, DISCARD, local_ele, temp, ax, no_norm, Summ, rhs,
						local_sino, osem_apu, N2, indices, elements, v_indices, idx, local_ind, N22);
					if (!DISCARD && !PRECOMP) {
						store_elements[tid + ind] = local_ele;
						store_indices[tid + ind] = local_ind;
						ind++;
					}
				}
				uu1 = xx;
				xx = 0;
				incr = 0;
				prev_local = 1.f;
				for (xx = alku_x2; xx >= 0; xx--) {
					double local_ele = compute_element_volume_3D(detectors, x_diff, y_diff, z_diff, kerroin, x_center[xx], y_center[yy2], z_center[zz]);
					if (local_ele >= bmax && incr > 0 && prev_local < local_ele) {
						if (xx == alku_x2 - 1) {
							breikki2 = true;
						}
						break;
					}
					else if (local_ele >= bmax) {
						incr++;
						prev_local = local_ele;
						continue;
					}
					incr = 1;
					prev_local = local_ele;
					if (local_ele < bmin)
						local_ele = Vmax;
					else
						local_ele = V[static_cast<size_t>(std::round((local_ele - bmin) * CC))];
					const uint32_t local_ind = compute_ind_orth_mfree_3D(static_cast<uint32_t>(xx), yy2 * N1, static_cast<uint32_t>(zz), NN, Nyx);
					computeIndices(RHS, SUMMA, OMP, PRECOMP, DISCARD, local_ele, temp, ax, no_norm, Summ, rhs,
						local_sino, osem_apu, N2, indices, elements, v_indices, idx, local_ind, N22);
					if (!DISCARD && !PRECOMP) {
						store_elements[tid + ind] = local_ele;
						store_indices[tid + ind] = local_ind;
						ind++;
					}
				}
				uu2 = xx;
				if (iu > 0) {
					if (ju > 0) {
						alku_x2 = uu2 + 1;
						alku_x1 = uu2 + 2;
					}
					else {
						alku_x1 = uu1 - 1;
						alku_x2 = uu1 - 2;
					}
				}
				else {
					if (ju > 0) {
						alku_x1 = uu1 - 1;
						alku_x2 = uu1 - 2;
					}
					else {
						alku_x2 = uu2 + 1;
						alku_x1 = uu2 + 2;
					}
				}
				if (breikki1 && breikki2) {
					if (yy2 == alku_y2)
						breikki4 = true;
					break;
				}
			}
			if (breikki3 && breikki4) {
				break;
			}
		}
		for (int zz = tempk - 1; zz >= loppu; zz--) {
			yy1 = 0;
			uu1 = 0;
			yy2 = 0;
			uu2 = 0;
			alku_y1 = tempj;
			alku_x1 = tempi;
			alku_y2 = tempj - 1;
			alku_x2 = tempi - 1;
			breikki1 = false;
			breikki2 = false;
			breikki3 = false;
			breikki4 = false;
			for (yy1 = alku_y1; yy1 < Ny; yy1++) {
				int xx = 0;
				int incr = 0;
				float prev_local = 1.f;
				for (xx = alku_x1; xx < Nx; xx++) {
					double local_ele = compute_element_volume_3D(detectors, x_diff, y_diff, z_diff, kerroin, x_center[xx], y_center[yy1], z_center[zz]);
					if (local_ele >= bmax && incr > 0 && prev_local < local_ele) {
						if (xx == alku_x1 + 1) {
							breikki1 = true;
						}
						break;
					}
					else if (local_ele >= bmax) {
						incr++;
						prev_local = local_ele;
						continue;
					}
					incr = 1;
					prev_local = local_ele;
					if (local_ele < bmin)
						local_ele = Vmax;
					else
						local_ele = V[static_cast<size_t>(std::round((local_ele - bmin) * CC))];
					const uint32_t local_ind = compute_ind_orth_mfree_3D(static_cast<uint32_t>(xx), yy1 * N1, static_cast<uint32_t>(zz), NN, Nyx);
					computeIndices(RHS, SUMMA, OMP, PRECOMP, DISCARD, local_ele, temp, ax, no_norm, Summ, rhs,
						local_sino, osem_apu, N2, indices, elements, v_indices, idx, local_ind, N22);
					if (!DISCARD && !PRECOMP) {
						store_elements[tid + ind] = local_ele;
						store_indices[tid + ind] = local_ind;
						ind++;
					}
				}
				uu1 = xx;
				xx = 0;
				incr = 0;
				prev_local = 1.f;
				for (xx = alku_x2; xx >= 0; xx--) {
					double local_ele = compute_element_volume_3D(detectors, x_diff, y_diff, z_diff, kerroin, x_center[xx], y_center[yy1], z_center[zz]);
					if (local_ele >= bmax && incr > 0 && prev_local < local_ele) {
						if (xx == alku_x2 - 1) {
							breikki2 = true;
						}
						break;
					}
					else if (local_ele >= bmax) {
						incr++;
						prev_local = local_ele;
						continue;
					}
					incr = 1;
					prev_local = local_ele;
					if (local_ele < bmin)
						local_ele = Vmax;
					else
						local_ele = V[static_cast<size_t>(std::round((local_ele - bmin) * CC))];
					const uint32_t local_ind = compute_ind_orth_mfree_3D(static_cast<uint32_t>(xx), yy1 * N1, static_cast<uint32_t>(zz), NN, Nyx);
					computeIndices(RHS, SUMMA, OMP, PRECOMP, DISCARD, local_ele, temp, ax, no_norm, Summ, rhs,
						local_sino, osem_apu, N2, indices, elements, v_indices, idx, local_ind, N22);
					if (!DISCARD && !PRECOMP) {
						store_elements[tid + ind] = local_ele;
						store_indices[tid + ind] = local_ind;
						ind++;
					}
				}
				uu2 = xx;
				if (iu > 0) {
					if (ju > 0) {
						alku_x1 = uu1 - 1;
						alku_x2 = uu1 - 2;
					}
					else {
						alku_x2 = uu2 + 1;
						alku_x1 = uu2 + 2;
					}
				}
				else {
					if (ju > 0) {
						alku_x2 = uu2 + 1;
						alku_x1 = uu2 + 2;
					}
					else {
						alku_x1 = uu1 - 1;
						alku_x2 = uu1 - 2;
					}
				}
				if (breikki1 && breikki2) {
					if (yy1 == alku_y1)
						breikki3 = true;
					break;
				}
			}
			breikki1 = false;
			breikki2 = false;
			alku_x1 = tempi;
			alku_x2 = tempi - 1;
			for (yy2 = alku_y2; yy2 >= 0; yy2--) {
				int xx = 0;
				int incr = 0;
				float prev_local = 1.f;
				for (xx = alku_x1; xx < Nx; xx++) {
					double local_ele = compute_element_volume_3D(detectors, x_diff, y_diff, z_diff, kerroin, x_center[xx], y_center[yy2], z_center[zz]);
					if (local_ele >= bmax && incr > 0 && prev_local < local_ele) {
						if (xx == alku_x1 + 1) {
							breikki1 = true;
						}
						break;
					}
					else if (local_ele >= bmax) {
						incr++;
						prev_local = local_ele;
						continue;
					}
					incr = 1;
					prev_local = local_ele;
					if (local_ele < bmin)
						local_ele = Vmax;
					else
						local_ele = V[static_cast<size_t>(std::round((local_ele - bmin) * CC))];
					const uint32_t local_ind = compute_ind_orth_mfree_3D(static_cast<uint32_t>(xx), yy2 * N1, static_cast<uint32_t>(zz), NN, Nyx);
					computeIndices(RHS, SUMMA, OMP, PRECOMP, DISCARD, local_ele, temp, ax, no_norm, Summ, rhs,
						local_sino, osem_apu, N2, indices, elements, v_indices, idx, local_ind, N22);
					if (!DISCARD && !PRECOMP) {
						store_elements[tid + ind] = local_ele;
						store_indices[tid + ind] = local_ind;
						ind++;
					}
				}
				uu1 = xx;
				xx = 0;
				incr = 0;
				prev_local = 1.f;
				for (xx = alku_x2; xx >= 0; xx--) {
					double local_ele = compute_element_volume_3D(detectors, x_diff, y_diff, z_diff, kerroin, x_center[xx], y_center[yy2], z_center[zz]);
					if (local_ele >= bmax && incr > 0 && prev_local < local_ele) {
						if (xx == alku_x2 - 1) {
							breikki2 = true;
						}
						break;
					}
					else if (local_ele >= bmax) {
						incr++;
						prev_local = local_ele;
						continue;
					}
					incr = 1;
					prev_local = local_ele;
					if (local_ele < bmin)
						local_ele = Vmax;
					else
						local_ele = V[static_cast<size_t>(std::round((local_ele - bmin) * CC))];
					const uint32_t local_ind = compute_ind_orth_mfree_3D(static_cast<uint32_t>(xx), yy2 * N1, static_cast<uint32_t>(zz), NN, Nyx);
					computeIndices(RHS, SUMMA, OMP, PRECOMP, DISCARD, local_ele, temp, ax, no_norm, Summ, rhs,
						local_sino, osem_apu, N2, indices, elements, v_indices, idx, local_ind, N22);
					if (!DISCARD && !PRECOMP) {
						store_elements[tid + ind] = local_ele;
						store_indices[tid + ind] = local_ind;
						ind++;
					}
				}
				uu2 = xx;
				if (iu > 0) {
					if (ju > 0) {
						alku_x2 = uu2 + 1;
						alku_x1 = uu2 + 2;
					}
					else {
						alku_x1 = uu1 - 1;
						alku_x2 = uu1 - 2;
					}
				}
				else {
					if (ju > 0) {
						alku_x1 = uu1 - 1;
						alku_x2 = uu1 - 2;
					}
					else {
						alku_x2 = uu2 + 1;
						alku_x1 = uu2 + 2;
					}
				}
				if (breikki1 && breikki2) {
					if (yy2 == alku_y2)
						breikki4 = true;
					break;
				}
			}
			if (breikki3 && breikki4) {
				break;
			}
		}

	}
}


// Calculate the denominator (forward projection) in the perpendicular case in orthogonal ray tracer (3D case)
void volume_distance_denominator_perpendicular_mfree_3D(const double* center1, const double center2, const double* z_center, double& temp,
	const bool d_attenuation_correction, const bool normalization, double& ax, const double d_b, const double d, const double d_d1, const uint32_t d_N1,
	const uint32_t d_N2, const uint32_t z_loop, const double* d_atten, const double* norm_coef, const double local_sino, const uint32_t d_N, const uint32_t d_NN, 
	const double* d_OSEM, Det detectors, const double xl, const double yl, const double zl, const double crystal_size_z, const uint32_t Nyx, 
	const uint32_t Nz, std::vector<double>& store_elements, std::vector<uint32_t>& store_indices, const uint32_t tid, uint32_t& ind, 
	double* elements, mwIndex* indices, const size_t lo, const bool PRECOMPUTE, const double global_factor, const double bmax, const double bmin, const double Vmax, 
	const double* V, const bool scatter, const double* scatter_coef, uint64_t N2) {

	//const uint32_t zz = z_loop * d_N2 * d_N1;
	const uint32_t apu = perpendicular_start(d_b, d, d_d1, d_N1);
	double jelppi = 0.;
	uint32_t hpk = N2;
	for (int32_t zz = static_cast<int32_t>(z_loop); zz >= 0; zz--) {
		for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
			double local_ele = compute_element_volume_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele >= bmax)
				break;
			if (local_ele < bmin)
				local_ele = Vmax;
			else
				local_ele = V[static_cast<size_t>(std::round((local_ele - bmin) * CC))];
			temp += (local_ele * d_N2);
			uint32_t local_ind = uu * d_N + zz * Nyx;
			if (!PRECOMPUTE) {
				store_indices[tid + ind] = local_ind;
				store_elements[tid + ind] = local_ele;
				ind++;
			}
			for (uint32_t kk = 0u; kk < d_N2; kk++) {
				if (d_attenuation_correction && uu == static_cast<int32_t>(apu) && zz == static_cast<int32_t>(z_loop))
					jelppi += (d_d1 * -d_atten[local_ind]);
				if (local_sino > 0.) {
					denominator_mfree(local_ele, ax, d_OSEM[local_ind]);
				}
				local_ind += d_NN;
			}
		}
		for (uint32_t uu = apu + 1; uu < d_N1; uu++) {
			double local_ele = compute_element_volume_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele >= bmax)
				break;
			if (local_ele < bmin)
				local_ele = Vmax;
			else
				local_ele = V[static_cast<size_t>(std::round((local_ele - bmin) * CC))];
			temp += (local_ele * d_N2);
			uint32_t local_ind = uu * d_N + zz * Nyx;
			if (!PRECOMPUTE) {
				store_indices[tid + ind] = local_ind;
				store_elements[tid + ind] = local_ele;
				ind++;
			}
			if (local_sino > 0.) {
				for (uint32_t kk = 0u; kk < d_N2; kk++) {
					denominator_mfree(local_ele, ax, d_OSEM[local_ind]);
					local_ind += d_NN;
				}
			}
		}
	}
	for (uint32_t zz = z_loop + 1u; zz < Nz; zz++) {
		for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
			double local_ele = compute_element_volume_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele >= bmax)
				break;
			temp += (local_ele * d_N2);
			uint32_t local_ind = uu * d_N + zz * Nyx;
			if (!PRECOMPUTE) {
				store_indices[tid + ind] = local_ind;
				store_elements[tid + ind] = local_ele;
				ind++;
			}
			if (local_sino > 0.) {
				for (uint32_t kk = 0u; kk < d_N2; kk++) {
					denominator_mfree(local_ele, ax, d_OSEM[local_ind]);
					local_ind += d_NN;
				}
			}
		}
		for (uint32_t uu = apu + 1; uu < d_N1; uu++) {
			double local_ele = compute_element_volume_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele >= bmax)
				break;
			temp += (local_ele * d_N2);
			uint32_t local_ind = uu * d_N + zz * Nyx;
			if (!PRECOMPUTE) {
				store_indices[tid + ind] = local_ind;
				store_elements[tid + ind] = local_ele;
				ind++;
			}
			if (local_sino > 0.) {
				for (uint32_t kk = 0u; kk < d_N2; kk++) {
					denominator_mfree(local_ele, ax, d_OSEM[local_ind]);
					local_ind += d_NN;
				}
			}
		}
	}
	temp = 1. / temp;
	if (d_attenuation_correction)
		temp *= exp(jelppi);
	if (normalization)
		temp *= norm_coef[lo];
	if (scatter)
		temp *= scatter_coef[lo];
	temp *= global_factor;
	if (PRECOMPUTE) {
		for (int32_t zz = static_cast<int32_t>(z_loop); zz >= 0; zz--) {
			for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
				double local_ele = compute_element_volume_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
				if (local_ele >= bmax)
					break;
				local_ele *= temp;
				uint32_t local_ind = uu * d_N + zz * Nyx;
				for (uint32_t kk = 0u; kk < d_N2; kk++) {
					indices[hpk] = local_ind;
					elements[hpk] = local_ele;
					local_ind += d_NN;
					hpk++;
				}
			}
			for (uint32_t uu = apu + 1; uu < d_N1; uu++) {
				double local_ele = compute_element_volume_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
				if (local_ele >= bmax)
					break;
				local_ele *= temp;
				uint32_t local_ind = uu * d_N + zz * Nyx;
				for (uint32_t kk = 0u; kk < d_N2; kk++) {
					indices[hpk] = local_ind;
					elements[hpk] = local_ele;
					local_ind += d_NN;
					hpk++;
				}
			}
		}
		for (uint32_t zz = z_loop + 1u; zz < Nz; zz++) {
			for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
				double local_ele = compute_element_volume_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
				if (local_ele >= bmax)
					break;
				local_ele *= temp;
				uint32_t local_ind = uu * d_N + zz * Nyx;
				for (uint32_t kk = 0u; kk < d_N2; kk++) {
					indices[hpk] = local_ind;
					elements[hpk] = local_ele;
					local_ind += d_NN;
					hpk++;
				}
			}
			for (uint32_t uu = apu + 1; uu < d_N1; uu++) {
				double local_ele = compute_element_volume_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
				if (local_ele >= bmax)
					break;
				local_ele *= temp;
				uint32_t local_ind = uu * d_N + zz * Nyx;
				for (uint32_t kk = 0u; kk < d_N2; kk++) {
					indices[hpk] = local_ind;
					elements[hpk] = local_ele;
					local_ind += d_NN;
					hpk++;
				}
			}
		}
	}
}

void volume_perpendicular_precompute(const uint32_t N1, const uint32_t N2, const uint32_t Nz, const double dd, const std::vector<double> vec,
	const double* center1, const double center2, const double* z_center, const double crystal_size_z, size_t& temp_koko, const Det detectors,
	const double xl, const double yl, const double zl, const uint32_t z_loop, const double bmax, const double bmin, const double Vmax,
	const double* V) {
	uint32_t apu = 0u;
	// Find the closest y-index value by finding the smallest y-distance between detector 2 and all the y-pixel coordinates
	for (size_t ii = 0ULL; ii < static_cast<size_t>(N2); ii++) {
		double temp = (vec[ii + 1ULL] - dd);
		if (temp > 0.) {
			apu = static_cast<uint32_t>(ii);
			break;
		}
	}
	uint32_t koko1 = 0u;
	uint32_t koko2 = 0u;
	for (int32_t zz = static_cast<int32_t>(z_loop); zz >= 0; zz--) {
		for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
			double local_ele = compute_element_volume_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele >= bmax)
				break;
			koko1++;
		}
		temp_koko += (koko1 * N2);
		for (uint32_t uu = apu + 1u; uu < N1; uu++) {
			double local_ele = compute_element_volume_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele >= bmax)
				break;
			koko2++;
		}
		temp_koko += (koko2 * N2);
		koko1 = 0u;
		koko2 = 0u;
	}
	for (uint32_t zz = z_loop + 1u; zz < Nz; zz++) {
		for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
			double local_ele = compute_element_volume_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele >= bmax)
				break;
			koko1++;
		}
		temp_koko += (koko1 * N2);
		for (uint32_t uu = apu + 1u; uu < N1; uu++) {
			double local_ele = compute_element_volume_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele >= bmax)
				break;
			koko2++;
		}
		temp_koko += (koko2 * N2);
		koko1 = 0u;
		koko2 = 0u;
	}
}