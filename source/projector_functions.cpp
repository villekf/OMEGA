/**************************************************************************
* Includes various functions required by the original Siddon, improved 
* Siddon and orthogonal distance based ray tracers.
* 
* References: 
* [1] Jacobs, F., Sundermann, E., De Sutter, B., Christiaens, M. Lemahieu, 
* I. (1998). A Fast Algorithm to Calculate the Exact Radiological Path 
* through a Pixel or Voxel Space. Journal of computing and information 
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
constexpr auto THR = 0.001;
constexpr auto H_THR = 1. - THR;

// Compute the orthogonal distance between a point and a line (in 3D)
double compute_element_orth_3D(Det detectors, const double xl, const double yl, const double zl, const double crystal_size_z, 
	const double xp, const double yp, const double zp) {

	double x1, y1, z1, x0, y0, z0;

	x0 = xp - detectors.xs;
	y0 = yp - detectors.ys;
	z0 = zp - detectors.zs;

	// Cross product
	x1 = yl * z0 - zl * y0;
	y1 = zl * x0 - xl * z0;
	z1 = xl * y0 - yl * x0;

	// Normalize the distance
	return (1. - norm(x1, y1, z1) / crystal_size_z);
}

// Compute the Euclidean norm of a vector
double norm(const double x, const double y, const double z) {
	return sqrt(x * x + y * y + z * z);
}

// Get the detector coordinates for the current raw list-mode measurement
void get_detector_coordinates_raw(const uint32_t det_per_ring, const double* x, const double* y, const double* z, Det& detectors,
	const uint16_t* L, const int ll, const uint32_t* pseudos, const uint32_t pRows) {
	uint32_t ps;
		const uint32_t detektorit1 = static_cast<uint32_t>(L[ll * 2]) - 1u;
		const uint32_t detektorit2 = static_cast<uint32_t>(L[ll * 2 + 1]) - 1u;

		const uint32_t loop1 = ((detektorit1) / det_per_ring);
		const uint32_t loop2 = ((detektorit2) / det_per_ring);

		if (loop1 == loop2) {
			detectors.zs = z[loop1];
			detectors.zd = detectors.zs;
		}
		else {
			detectors.zs = z[loop1];
			detectors.zd = z[loop2];
		}

		if (loop1 >= pseudos[0]) {
			ps = 1u;
			for (uint32_t kk = 0u; kk < pRows; kk++) {
				if (kk + 1u < pRows) {
					if (loop1 >= pseudos[kk] && loop1 < pseudos[kk + 1]) {
						detectors.zs = z[loop1 + ps];
						break;
					}
					else
						ps++;
				}
				else {
					if (loop1 >= pseudos[kk])
						detectors.zs = z[loop1 + ps];
				}
			}
		}

		if (loop2 >= pseudos[0]) {
			ps = 1u;
			for (uint32_t kk = 0u; kk < pRows; kk++) {
				if (kk + 1 < pRows) {
					if (loop2 >= pseudos[kk] && loop2 < pseudos[kk + 1]) {
						detectors.zd = z[loop2 + ps];
						break;
					}
					else
						ps++;
				}
				else {
					if (loop2 >= pseudos[kk])
						detectors.zd = z[loop2 + ps];
				}
			}
		}

		detectors.xs = x[detektorit1 - det_per_ring * (loop1)];
		detectors.xd = x[detektorit2 - det_per_ring * (loop2)];
		detectors.ys = y[detektorit1 - det_per_ring * (loop1)];
		detectors.yd = y[detektorit2 - det_per_ring * (loop2)];
	
}

// Get the detector coordinates for the current raw list-mode measurement (multi-ray)
void get_detector_coordinates_raw_N(const uint32_t det_per_ring, const double* x, const double* y, const double* z, Det& detectors,
	const uint16_t* L, const int ll, const uint32_t* pseudos, const uint32_t pRows, const uint16_t lor, const double cr_pz) {
	uint32_t ps;
	const uint32_t detektorit1 = static_cast<uint32_t>(L[ll * 2]) - 1u;
	const uint32_t detektorit2 = static_cast<uint32_t>(L[ll * 2 + 1]) - 1u;

	const uint32_t loop1 = ((detektorit1) / det_per_ring);
	const uint32_t loop2 = ((detektorit2) / det_per_ring);

	if (loop1 == loop2) {
		detectors.zs = z[loop1];
		detectors.zd = detectors.zs;
	}
	else {
		detectors.zs = z[loop1];
		detectors.zd = z[loop2];
	}

	if (loop1 >= pseudos[0]) {
		ps = 1u;
		for (uint32_t kk = 0u; kk < pRows; kk++) {
			if (kk + 1u < pRows) {
				if (loop1 >= pseudos[kk] && loop1 < pseudos[kk + 1]) {
					detectors.zs = z[loop1 + ps];
					break;
				}
				else
					ps++;
			}
			else {
				if (loop1 >= pseudos[kk])
					detectors.zs = z[loop1 + ps];
			}
		}
	}

	if (loop2 >= pseudos[0]) {
		ps = 1u;
		for (uint32_t kk = 0u; kk < pRows; kk++) {
			if (kk + 1 < pRows) {
				if (loop2 >= pseudos[kk] && loop2 < pseudos[kk + 1]) {
					detectors.zd = z[loop2 + ps];
					break;
				}
				else
					ps++;
			}
			else {
				if (loop2 >= pseudos[kk])
					detectors.zd = z[loop2 + ps];
			}
		}
	}

	if (lor == 1) {
		detectors.xs = x[detektorit1 - det_per_ring * (loop1)];
		detectors.xd = x[detektorit2 - det_per_ring * (loop2)];
		detectors.ys = y[detektorit1 - det_per_ring * (loop1)];
		detectors.yd = y[detektorit2 - det_per_ring * (loop2)];
	}
	else if (lor == 3u || lor == 5u) {
		detectors.xs = x[detektorit1 - det_per_ring * (loop1) + det_per_ring * 2u];
		detectors.xd = x[detektorit2 - det_per_ring * (loop2) + det_per_ring * 2u];
		detectors.ys = y[detektorit1 - det_per_ring * (loop1) + det_per_ring * 2u];
		detectors.yd = y[detektorit2 - det_per_ring * (loop2) + det_per_ring * 2u];
		if (lor == 3u) {
			detectors.zs -= cr_pz;
			detectors.zd -= cr_pz;
		}
		else {
			detectors.zs += cr_pz;
			detectors.zd += cr_pz;
		}
	}
	else {
		detectors.xs = x[detektorit1 - det_per_ring * (loop1) + det_per_ring];
		detectors.xd = x[detektorit2 - det_per_ring * (loop2) + det_per_ring];
		detectors.ys = y[detektorit1 - det_per_ring * (loop1) + det_per_ring];
		detectors.yd = y[detektorit2 - det_per_ring * (loop2) + det_per_ring];
		if (lor == 2u) {
			detectors.zs -= cr_pz;
			detectors.zd -= cr_pz;
		}
		else {
			detectors.zs += cr_pz;
			detectors.zd += cr_pz;
		}
	}

}

// Get the detector coordinates for the current sinogram bin, no precomputations performed beforehand
void get_detector_coordinates_noalloc(const double* x, const double* y, const double* z, const uint32_t size_x, Det& detectors, int& ll, 
	const uint32_t* index, int& lz, const uint32_t TotSinos, uint32_t oo) {
	// Sinogram data

	ll = (index[oo] - 1u) % size_x;
	lz = (index[oo] - 1u) / size_x;

	detectors.xs = x[ll];
	detectors.xd = x[ll + size_x];
	detectors.ys = y[ll];
	detectors.yd = y[ll + size_x];
	detectors.zs = z[lz];
	detectors.zd = z[lz + TotSinos];
}

// Get the detector coordinates for the current sinogram bin, precomputed data
void get_detector_coordinates(const double* x, const double* y, const double* z, const uint32_t size_x, Det& detectors, const uint32_t* xy_index,
	const uint16_t* z_index, const uint32_t TotSinos, const uint32_t oo) {
	// Sinogram data
	if (xy_index[oo] >= size_x) {
		detectors.xs = x[xy_index[oo]];
		detectors.xd = x[xy_index[oo] - size_x];
		detectors.ys = y[xy_index[oo]];
		detectors.yd = y[xy_index[oo] - size_x];
	}
	else {
		detectors.xs = x[xy_index[oo]];
		detectors.xd = x[xy_index[oo] + size_x];
		detectors.ys = y[xy_index[oo]];
		detectors.yd = y[xy_index[oo] + size_x];
	}
	if (z_index[oo] >= TotSinos) {
		detectors.zs = z[z_index[oo]];
		detectors.zd = z[z_index[oo] - TotSinos];
	}
	else {
		detectors.zs = z[z_index[oo]];
		detectors.zd = z[z_index[oo] + TotSinos];
	}
}

// Get the detector coordinates for the current sinogram bin, precomputed data
void get_detector_coordinates_mr(const double* x, const double* y, const double* z, const uint32_t size_x, Det& detectors, const uint32_t* xy_index,
	const uint16_t* z_index, const uint32_t TotSinos, const uint32_t oo, const uint16_t lor, const double cr_pz) {
	// Sinogram data
	if (lor == 1u) {
		if (xy_index[oo] >= size_x) {
			detectors.xs = x[xy_index[oo]];
			detectors.xd = x[xy_index[oo] - size_x];
			detectors.ys = y[xy_index[oo]];
			detectors.yd = y[xy_index[oo] - size_x];
		}
		else {
			detectors.xs = x[xy_index[oo]];
			detectors.xd = x[xy_index[oo] + size_x];
			detectors.ys = y[xy_index[oo]];
			detectors.yd = y[xy_index[oo] + size_x];
		}
		if (z_index[oo] >= TotSinos) {
			detectors.zs = z[z_index[oo]];
			detectors.zd = z[z_index[oo] - TotSinos];
		}
		else {
			detectors.zs = z[z_index[oo]];
			detectors.zd = z[z_index[oo] + TotSinos];
		}
	}
	else if (lor == 3u || lor == 5u) {
		if (xy_index[oo] >= size_x) {
			detectors.xs = x[xy_index[oo] + size_x * 4u];
			detectors.xd = x[xy_index[oo] + size_x * 3u];
			detectors.ys = y[xy_index[oo] + size_x * 4u];
			detectors.yd = y[xy_index[oo] + size_x * 3u];
		}
		else {
			detectors.xs = x[xy_index[oo] + size_x * 4u];
			detectors.xd = x[xy_index[oo] + size_x * 5u];
			detectors.ys = y[xy_index[oo] + size_x * 4u];
			detectors.yd = y[xy_index[oo] + size_x * 5u];
		}
		if (lor == 3u) {
			if (z_index[oo] >= TotSinos) {
				detectors.zs = z[z_index[oo]] - cr_pz;
				detectors.zd = z[z_index[oo] - TotSinos] - cr_pz;
			}
			else {
				detectors.zs = z[z_index[oo]] - cr_pz;
				detectors.zd = z[z_index[oo] + TotSinos] - cr_pz;
			}
		}
		else {
			if (z_index[oo] >= TotSinos) {
				detectors.zs = z[z_index[oo]] + cr_pz;
				detectors.zd = z[z_index[oo] - TotSinos] + cr_pz;
			}
			else {
				detectors.zs = z[z_index[oo]] + cr_pz;
				detectors.zd = z[z_index[oo] + TotSinos] + cr_pz;
			}
		}
	}
	else {
		if (xy_index[oo] >= size_x) {
			detectors.xs = x[xy_index[oo] + size_x * 2u];
			detectors.xd = x[xy_index[oo] + size_x];
			detectors.ys = y[xy_index[oo] + size_x * 2u];
			detectors.yd = y[xy_index[oo] + size_x];
		}
		else {
			detectors.xs = x[xy_index[oo] + size_x * 2u];
			detectors.xd = x[xy_index[oo] + size_x * 3u];
			detectors.ys = y[xy_index[oo] + size_x * 2u];
			detectors.yd = y[xy_index[oo] + size_x * 3u];
		}
		if (lor == 2u) {
			if (z_index[oo] >= TotSinos) {
				detectors.zs = z[z_index[oo]] - cr_pz;
				detectors.zd = z[z_index[oo] - TotSinos] - cr_pz;
			}
			else {
				detectors.zs = z[z_index[oo]] - cr_pz;
				detectors.zd = z[z_index[oo] + TotSinos] - cr_pz;
			}
		}
		else {
			if (z_index[oo] >= TotSinos) {
				detectors.zs = z[z_index[oo]] + cr_pz;
				detectors.zd = z[z_index[oo] - TotSinos] + cr_pz;
			}
			else {
				detectors.zs = z[z_index[oo]] + cr_pz;
				detectors.zd = z[z_index[oo] + TotSinos] + cr_pz;
			}
		}
	}
}

// Get the current ring in axial direction
uint32_t z_ring(const double zmax, const double zs, const double NSlices) {
	return static_cast<uint32_t>((zs / zmax)*(NSlices - 1.));
}

// Compute the initial and maximum voxel indices and the direction of the ray, detector greater than source
void d_g_s(const double tmin, const double t_min, const double tmax, const double t_max, uint32_t& v_min, uint32_t& v_max, double& t_0, int32_t& v_u, 
	const double diff, const double b, const double d, const double s, const uint32_t N) {

	if (tmin == t_min)
		// (11)
		v_min = 1u;
	else {
		// (2) and (19)
		const double p_t = s + tmin * (diff);
		// (12)
		v_min = static_cast<uint32_t>(ceil((p_t - b) / d));
	}
	//if (tmax == t_max)
	//	// (13)
	//	v_max = N;
	//else {
	//	// (2) and (19)
	//	p_t = s + tmax * (diff);
	//	// (14)
	//	v_max = static_cast<int>(floor((p_t - b) / d));
	//}
	// (9)
	//tx0 = (bx + static_cast<double>(imin) * dx - xs) / (x_diff);
	t_0 += static_cast<double>(v_min) * d / (diff);
	//  (29)
	v_u = 1;
}

// Compute the initial and maximum voxel indices and the direction of the ray, source greater than detector
void s_g_d(const double tmin, const double t_min, const double tmax, const double t_max, uint32_t& v_min, uint32_t& v_max, double& t_0, int32_t& v_u, 
	const double diff, const double b, const double d, const double s, const uint32_t N) {

	if (tmin == t_min)
		// (15)
		v_max = N - 1u;
	else {
		// (2) and (19)
		const double p_t = s + tmin * (diff);
		// (16)
		v_max = static_cast<uint32_t>(floor((p_t - b) / d));
	}
	//if (tmax == t_max)
	//	// (17)
	//	v_min = 0;
	//else {
	//	// (2) and (19)
	//	p_t = s + tmax * (diff);
	//	// (18)
	//	v_min = static_cast<int>(ceil((p_t - b) / d));
	//}
	// (9)
	//tx0 = (bx + static_cast<double>(imax) * dx - xs) / (x_diff);
	t_0 += static_cast<double>(v_max) * d / (diff);
	// (29)
	v_u = -1;
}

// Compute the initial and maximum voxel indices and the direction of the ray, detector greater than source
void d_g_s_precomp(const double tmin, const double t_min, const double tmax, const double t_max, uint32_t& v_min, uint32_t& v_max, double& t_0, 
	int32_t& v_u, const double diff, const double b, const double d, const double s, const uint32_t N) {

	if (tmin == t_min)
		// (11)
		v_min = 1u;
	else {
		// (2) and (19)
		const double p_t = s + tmin * (diff);
		// (12)
		v_min = static_cast<uint32_t>(ceil((p_t - b) / d));
	}
	if (tmax == t_max)
		// (13)
		v_max = N;
	else {
		// (2) and (19)
		const double p_t = s + tmax * (diff);
		// (14)
		v_max = static_cast<uint32_t>(floor((p_t - b) / d));
	}
	// (9)
	//tx0 = (bx + static_cast<double>(imin) * dx - xs) / (x_diff);
	t_0 += static_cast<double>(v_min) * d / (diff);
	//  (29)
	v_u = 1;
}

// Compute the initial and maximum voxel indices and the direction of the ray, source greater than detector
void s_g_d_precomp(const double tmin, const double t_min, const double tmax, const double t_max, uint32_t& v_min, uint32_t& v_max, double& t_0, 
	int32_t& v_u, const double diff, const double b, const double d, const double s, const uint32_t N) {

	if (tmin == t_min)
		// (15)
		v_max = N - 1u;
	else {
		// (2) and (19)
		const double p_t = s + tmin * (diff);
		// (16)
		v_max = static_cast<uint32_t>(floor((p_t - b) / d));
	}
	if (tmax == t_max)
		// (17)
		v_min = 0u;
	else {
		// (2) and (19)
		const double p_t = s + tmax * (diff);
		// (18)
		v_min = static_cast<uint32_t>(ceil((p_t - b) / d));
	}
	// (9)
	//tx0 = (bx + static_cast<double>(imax) * dx - xs) / (x_diff);
	t_0 += static_cast<double>(v_max) * d / (diff);
	// (29)
	v_u = -1;
}

// Compute the distance that the ray traverses in the current voxel
double pixel_value(const double t, const double tc, const double L) {
	return (t - tc) * L;
}

void compute_attenuation(double& tc, double& jelppi, const double LL, const double t0, const int tempi, const int tempj, const int tempk, 
	const uint32_t Nx, const uint32_t Nyx, const double* atten) {
	jelppi += (pixel_value(t0, tc, LL) * -atten[tempi + tempj * Nx + Nyx * tempk]);
	tc = t0;
}

// Correct for attenuation, vector data
void att_corr_vec(const std::vector<double> templ_ijk, const std::vector<uint32_t> temp_koko, const double* atten, double& temp, const size_t Np) {

	double jelppi = 0.;

	for (uint32_t iii = 0u; iii < Np; iii++) {
		jelppi += templ_ijk[iii] * -atten[temp_koko[iii]];
	}
	temp = std::exp(jelppi) * temp;
}

// Correct for attenuation, vector data, precomputed
void att_corr_vec_precomp(const double* elements, const double* atten, const mwIndex* indices, const size_t Np, const uint64_t N2, double& temp) {

	double jelppi = 0.;

	for (uint32_t iii = 0u; iii < Np; iii++) {
		jelppi += (elements[N2 + iii]) * -atten[indices[N2 + iii]];
	}
	temp = exp(jelppi) * temp;
}

// Correct for attenuation, scalar data
double att_corr_scalar(double templ_ijk, uint32_t tempk, const double* atten, double& temp, const uint32_t N1, const uint32_t N) {

	double jelppi = 0.;
	for (uint32_t iii = 0u; iii < N1; iii++) {
		jelppi += templ_ijk * -atten[tempk + iii * N];
	}
	temp *= std::exp(jelppi);
	return jelppi;
}

// Correct for attenuation, orthogonal distance based ray tracer
void att_corr_scalar_orth(uint32_t tempk, const double* atten, double& temp, const uint32_t N1, const uint32_t N2, const double d) {

	double jelppi = 0.;
	for (uint32_t ii = 0u; ii < N2; ii++) {
		jelppi += d * -atten[tempk + ii * N1];
	}
	temp *= std::exp(jelppi);
}

// Compute the probability for one emission in perpendicular detector case
double perpendicular_elements(const uint32_t N, const double dd, const std::vector<double> vec, const double d, const uint32_t z_ring, 
	const uint32_t N1, const uint32_t N2, const double* atten, const double* norm_coef, const bool attenuation_correction, const bool normalization, 
	uint32_t& tempk, const uint32_t NN, const uint32_t lo) {
	uint32_t apu = 0u;
	// Find the closest y-index value by finding the smallest y-distance between detector 2 and all the y-pixel coordinates
	for (size_t ii = 0ULL; ii < static_cast<size_t>(N2); ii++) {
		double temp = (vec[ii + 1ULL] - dd);
		if (temp > 0.) {
			apu = static_cast<uint32_t>(ii);
			break;
		}
	}
	tempk = apu * N + z_ring * N1 * N2;
	double temp = d * static_cast<double>(N2);
	// Probability
	temp = 1. / temp;

	// Correct for attenuation if applicable
	if (attenuation_correction)
		att_corr_scalar(d, tempk, atten, temp, N2, NN);
	if (normalization)
		temp *= norm_coef[lo];


	return d * temp;
}

// Compute the probability for one emission in perpendicular detector case (multi-ray)
double perpendicular_elements_multiray(const uint32_t N, const double dd, const std::vector<double> vec, const double d, const uint32_t z_ring,
	const uint32_t N1, const uint32_t N2, const double* atten, const bool attenuation_correction, int32_t& tempk, const uint32_t NN, double& jelppi) {
	uint32_t apu = 0u;
	// Find the closest y-index value by finding the smallest y-distance between detector 2 and all the y-pixel coordinates
	for (size_t ii = 0ULL; ii < static_cast<size_t>(N2); ii++) {
		double temp = (vec[ii + 1ULL] - dd);
		if (temp > 0.) {
			apu = static_cast<uint32_t>(ii);
			break;
		}
	}
	tempk = apu * N + z_ring * N1 * N2;
	double temp = 0.;

	// Correct for attenuation if applicable
	if (attenuation_correction)
		jelppi += att_corr_scalar(d, tempk, atten, temp, N2, NN);


	return d * static_cast<double>(N2);
}

// Compute the first voxel index
int32_t voxel_index(const double pt, const double diff, const double d, const double apu) {
	return static_cast<int32_t>(floor((pt * diff - apu) / d));
}

// compute the orthogonal distance, perpendicular detectors
void orth_perpendicular_precomputed(const uint32_t N, const double dd, const std::vector<double> vec, const double d, const uint32_t z_ring, 
	const uint32_t N1, const uint32_t N2, std::vector<uint32_t>& indices, std::vector<double>& elements, const double xd, const double xs, 
	const double yd, const double ys, const double diff2, const double* center1, const double crystal_size, const double length, 
	std::vector<double>& vec1, std::vector<double>& vec2, double& temp, uint32_t& tempk) {
	uint32_t apu = 0u;
	// Find the closest y-index value by finding the smallest y-distance between detector 2 and all the y-pixel coordinates
	for (size_t ii = 0ULL; ii < static_cast<size_t>(N2); ii++) {
		double temp = (vec[ii + 1ULL] - dd);
		if (temp > 0.) {
			apu = static_cast<uint32_t>(ii);
			break;
		}
	}
	tempk = apu * N + z_ring * N1 * N2;
	const double kerroin = xd * ys - yd * xs;
	for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
		const double d_ort = 1. - std::fabs(kerroin - diff2 * center1[uu]) / length;
		if (d_ort <= THR)
			break;
		vec1.emplace_back(d_ort);
		temp += d_ort;
	}
	for (uint32_t ky = apu + 1u; ky < N2; ky++) {
		const double d_ort = 1. - std::fabs(kerroin - diff2 * center1[ky]) / length;
		if (d_ort <= THR)
			break;
		vec2.emplace_back(d_ort);
		temp += d_ort;
	}
	temp += temp * static_cast<double>(N2 - 1);
	// Probability
	temp = 1. / temp;
}

// compute the orthogonal distance, perpendicular detectors
void orth_perpendicular(const uint32_t N, const double dd, const std::vector<double> vec, const double d, const uint32_t z_ring, const uint32_t N1, 
	const uint32_t N2, const double diff2, const double* center1, const double kerroin, const double length, std::vector<double>& vec1, 
	std::vector<double>& vec2, double& temp, uint32_t& tempk) {
	uint32_t apu = 0u;
	// Find the closest y-index value by finding the smallest y-distance between detector 2 and all the y-pixel coordinates
	for (size_t ii = 0ULL; ii < static_cast<size_t>(N2); ii++) {
		double temp = (vec[ii + 1ULL] - dd);
		if (temp > 0.) {
			apu = static_cast<uint32_t>(ii);
			break;
		}
	}
	tempk = apu * N + z_ring * N1 * N2;
	for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
		const double d_ort = compute_element_orth_mfree(kerroin, diff2, center1[uu], length);
		if (d_ort <= THR)
			break;
		vec1.emplace_back(d_ort);
		temp += d_ort;
	}
	for (uint32_t uu = apu + 1u; uu < N2; uu++) {
		const double d_ort = compute_element_orth_mfree(kerroin, diff2, center1[uu], length);
		if (d_ort <= THR)
			break;
		vec2.emplace_back(d_ort);
		temp += d_ort;
	}
	temp += temp * static_cast<double>(N2 - 1u);
	// Probability
	temp = 1. / temp;
}

// compute the orthogonal distance, perpendicular detectors
void orth_perpendicular_np_3D(const double dd, const std::vector<double> vec, const uint32_t z_ring, const uint32_t N1, const uint32_t N2, 
	const uint32_t Nz, const uint32_t Nyx, const uint32_t d_N, const uint32_t d_NN, const Det detectors, const double xl, const double yl, 
	const double zl, const double* center1, const double center2, const double* z_center, const double crystal_size_z, int& hpk, double& temp, 
	uint32_t& tempk, std::vector<uint32_T>& indices, std::vector<double>& elements) {
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
			double d_ort = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (d_ort <= THR)
				break;
			uint32_t local_ind = uu * d_N + zz * Nyx;
			for (uint32_t kk = 0u; kk < N1; kk++) {
				indices.emplace_back(local_ind + kk * d_NN);
				elements.emplace_back(d_ort);
				hpk++;
			}
			temp += d_ort * static_cast<double>(N2);
		}
		for (uint32_t uu = apu + 1; uu < N2; uu++) {
			double d_ort = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (d_ort <= THR)
				break;
			uint32_t local_ind = uu * d_N + zz * Nyx;
			temp += d_ort;
			for (uint32_t kk = 0u; kk < N1; kk++) {
				indices.emplace_back(local_ind + kk * d_NN);
				elements.emplace_back(d_ort);
				hpk++;
			}
			temp += d_ort * static_cast<double>(N2);
		}
	}
	for (uint32_t zz = z_ring + 1u; zz < Nz; zz++) {
		for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
			double d_ort = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (d_ort <= THR)
				break;
			uint32_t local_ind = uu * d_N + zz * Nyx;
			temp += d_ort;
			for (uint32_t kk = 0u; kk < N1; kk++) {
				indices.emplace_back(local_ind + kk * d_NN);
				elements.emplace_back(d_ort);
				hpk++;
			}
			temp += d_ort * static_cast<double>(N2);
		}
		for (uint32_t uu = apu + 1; uu < N2; uu++) {
			double d_ort = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (d_ort <= THR)
				break;
			uint32_t local_ind = uu * d_N + zz * Nyx;
			temp += d_ort;
			for (uint32_t kk = 0u; kk < N1; kk++) {
				indices.emplace_back(local_ind + kk * d_NN);
				elements.emplace_back(d_ort);
				hpk++;
			}
			temp += d_ort * static_cast<double>(N2);
		}
	}
	// Probability
	temp = 1. / temp;
}


// compute the orthogonal distance, perpendicular detectors
void orth_perpendicular_3D(const double dd, const std::vector<double> vec, const uint32_t z_ring, const uint32_t N1, const uint32_t N2, 
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
			double d_ort = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (d_ort <= THR)
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
			double d_ort = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (d_ort <= THR)
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
			double d_ort = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (d_ort <= THR)
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
			double d_ort = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (d_ort <= THR)
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

// compute the orthogonal distance, perpendicular detectors
void orth_perpendicular_precompute(const uint32_t N1, const uint32_t N2, const double dd, const std::vector<double> vec, const double diff2, 
	const double* center1, const double kerroin, const double length, uint16_t &temp_koko) {
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
	for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
		const double d_ort = std::fabs(kerroin + diff2 * center1[uu]) / length;
		if (d_ort >= H_THR)
			break;
		koko1++;
	}
	temp_koko += koko1 * N1;
	for (uint32_t uu = apu + 1u; uu < N2; uu++) {
		const double d_ort = std::fabs(kerroin + diff2 * center1[uu]) / length;
		if (d_ort >= H_THR)
			break;
		koko2++;
	}
	temp_koko += koko2 * N1;
}

// compute the orthogonal distance, perpendicular detectors
void orth_perpendicular_precompute_3D(const uint32_t N1, const uint32_t N2, const uint32_t Nz, const double dd, const std::vector<double> vec, 
	const double* center1, const double center2, const double* z_center, const double crystal_size_z, uint16_t& temp_koko, const Det detectors,
	const double xl, const double yl, const double zl, const uint32_t z_loop) {
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
			const double d_ort = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (d_ort <= THR)
				break;
			koko1++;
		}
		temp_koko += (koko1 * N2);
		for (uint32_t uu = apu + 1u; uu < N1; uu++) {
			const double d_ort = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (d_ort <= THR)
				break;
			koko2++;
		}
		temp_koko += (koko2 * N2);
		koko1 = 0u;
		koko2 = 0u;
	}
	for (uint32_t zz = z_loop + 1u; zz < Nz; zz++) {
		for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
			const double d_ort = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (d_ort <= THR)
				break;
			koko1++;
		}
		temp_koko += (koko1 * N2);
		for (uint32_t uu = apu + 1u; uu < N1; uu++) {
			const double d_ort = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (d_ort <= THR)
				break;
			koko2++;
		}
		temp_koko += (koko2 * N2);
		koko1 = 0u;
		koko2 = 0u;
	}
}

void orth_distance_precompute(const int32_t tempi, const uint32_t Nx, const double y_diff, const double x_diff, const double* x_center, 
	const double y_center, const double kerroin, const double length_, uint16_t& temp_koko) {

	const double diff = kerroin - x_diff * y_center;
	for (int32_t uu = tempi; uu >= 0; uu--) {
		const double d_ort = compute_element_orth_mfree(diff, y_diff, x_center[uu], length_);
		if (d_ort <= THR)
			break;
		temp_koko++;
	}
	for (uint32_t uu = static_cast<uint32_t>(tempi) + 1u; uu < Nx; uu++) {
		const double d_ort = compute_element_orth_mfree(diff, y_diff, x_center[uu], length_);
		if (d_ort <= THR)
			break;
		temp_koko++;
	}
}

void orth_distance_precompute_3D(const int32_t tempi, const uint32_t N, const uint32_t Nz, const double y_diff, const double x_diff, const double z_diff, 
	const double* x_center, const double y_center, const double* z_center, const double crystal_size_z, const Det detectors, 
	uint16_t& temp_koko, const int32_t tempk, const uint32_t lo, const int32_t n_tempk, const int32_t dec, const int32_t decx) {

	bool loppu1 = true;
	bool loppu2 = true;
	uint8_t loppu3 = 0u;
	uint8_t loppu4 = 0u;
	bool pass1 = true;
	bool pass = false;
	int32_t alku = tempk;
	int32_t alku2 = tempk;
	int32_t alku3 = tempk;
	for (int32_t uu = tempi; uu >= 0; uu--) {
		//for (int32_t zz = tempk; zz >= 0; zz--) {
		for (int32_t zz = alku; zz >= 0; zz--) {
			double local_ele = compute_element_orth_3D(detectors, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
			if (local_ele <= THR) {
				if (loppu1 && zz >= n_tempk - dec && !(pass && alku3 == tempk)) {
					continue;
				}
				else {
					if (loppu1 && loppu3 == 1u && loppu4 == 0u)
						pass1 = false;
					break;
				}
			}
			temp_koko++;
			if (loppu1) {
				alku = zz;
				loppu1 = false;
				loppu3 = 1u;
				if (zz == tempk || loppu4 == 1u)
					pass = true;
			}
			alku3 = zz;
		}
		if (pass) {
			for (int32_t zz = alku2 + 1; zz < Nz; zz++) {
				//for (uint32_t zz = static_cast<uint32_t>(tempk) + 1u; zz < Nz; zz++) {
				double local_ele = compute_element_orth_3D(detectors, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
				if (local_ele <= THR) {
					if (loppu2 && zz <= tempk + 1 + decx) {
						continue;
					}
					else {
						if (loppu2 && loppu4 == 1u) {
							loppu4 = 0u;
							pass = false;
						}
						break;
					}
				}
				temp_koko++;
				if (loppu2) {
					alku2 = zz - 1;
					loppu2 = false;
					loppu4 = 1u;
					pass1 = true;
				}
			}
		}
		loppu1 = true;
		loppu2 = true;
		if (pass1) {
			continue;
		}
		else {
			break;
		}
	}
	loppu3 = 0u;
	pass1 = true;
	pass = false;
	loppu4 = 0u;
	alku2 = tempk;
	alku = tempk;
	alku3 = tempk;
	for (uint32_t uu = static_cast<uint32_t>(tempi) + 1u; uu < N; uu++) {
		for (int32_t zz = alku; zz >= 0; zz--) {
			//for (int32_t zz = tempk; zz >= 0; zz--) {
			double local_ele = compute_element_orth_3D(detectors, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
			if (local_ele <= THR) {
				if (loppu1 && zz >= n_tempk - dec && !(pass && alku3 == tempk)) {
					continue;
				}
				else {
					if (loppu1 && loppu3 == 1u && loppu4 == 0u)
						pass1 = false;
					break;
				}
			}
			temp_koko++;
			if (loppu1) {
				alku = zz;
				loppu1 = false;
				loppu3 = 1u;
				if (zz == tempk || loppu4 == 1u)
					pass = true;
			}
			alku3 = zz;
		}
		if (pass) {
			for (int32_t zz = alku2 + 1; zz < Nz; zz++) {
				//for (uint32_t zz = static_cast<uint32_t>(tempk) + 1u; zz < Nz; zz++) {
				double local_ele = compute_element_orth_3D(detectors, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
				if (local_ele <= THR) {
					if (loppu2 && zz <= tempk + 1 + decx) {
						continue;
					}
					else {
						if (loppu2 && loppu4 == 1u) {
							loppu4 = 0u;
							pass = false;
						}
						break;
					}
				}
				temp_koko++;
				if (loppu2) {
					alku2 = zz - 1;
					loppu2 = false;
					loppu4 = 1u;
					pass1 = true;
				}
			}
		}
		loppu1 = true;
		loppu2 = true;
		if (pass1) {
			continue;
		}
		else {
			break;
		}
	}

}

// Compute whether the ray intercepts the FOV (TYPE == 0), compute the voxel indices of the first voxel intercepted, 
// total number of voxels traversed and coefficients needed for Siddon's algorithm (2D)
bool siddon_pre_loop_2D(const double b1, const double b2, const double diff1, const double diff2, const double max1, const double max2, 
	const double d1, const double d2, const uint32_t N1, const uint32_t N2, int32_t& temp1, int32_t& temp2, double& t1u, double& t2u, uint32_t& Np,
	const int TYPE, const double ys, const double xs, const double yd, const double xd, double& tc, int32_t& u1, int32_t& u2, double& t10, double& t20) {
	// If neither x- nor y-directions are perpendicular
// Correspond to the equations (9) and (10) from reference [1]
	const double apu_tx = b1 - xs;
	const double apu_ty = b2 - ys;
	t10 = (apu_tx) / (diff1);
	t20 = (apu_ty) / (diff2);
	const double txback = (max1 - xs) / (diff1);
	const double tyback = (max2 - ys) / (diff2);

	// Equations (5-8)
	const double txmin = std::min(t10, txback);
	const double txmax = std::max(t10, txback);
	const double tymin = std::min(t20, tyback);
	const double tymax = std::max(t20, tyback);

	// (3-4)
	tc = std::max(txmin, tymin);
	const double tmax = std::min(txmax, tymax);

	uint32_t imin, imax, jmin, jmax;

	if (TYPE == 0) {
	// If true, then the ray/LOR does not intersect the pixel space --> continue to the next LOR
		if (tc >= tmax) {
			return true;
		}

		// (11-14)
		if (xs < xd)
			d_g_s_precomp(tc, txmin, tmax, txmax, imin, imax, t10, u1, diff1, b1, d1, xs, N1);
		// (15-18)
		else
			s_g_d_precomp(tc, txmin, tmax, txmax, imin, imax, t10, u1, diff1, b1, d1, xs, N1);

		//Same as above
		if (ys < yd)
			d_g_s_precomp(tc, tymin, tmax, tymax, jmin, jmax, t20, u2, diff2, b2, d2, ys, N2);
		else
			s_g_d_precomp(tc, tymin, tmax, tymax, jmin, jmax, t20, u2, diff2, b2, d2, ys, N2);

		Np = imax + 1u + jmax + 1u - imin - jmin;

		//tc = tmin;
	}
	else {
		// (11-14)
		if (xs < xd)
			d_g_s(tc, txmin, tmax, txmax, imin, imax, t10, u1, diff1, b1, d1, xs, N1);
		// (15-18)
		else
			s_g_d(tc, txmin, tmax, txmax, imin, imax, t10, u1, diff1, b1, d1, xs, N1);
		//Same as above
		if (ys < yd)
			d_g_s(tc, tymin, tmax, tymax, jmin, jmax, t20, u2, diff2, b2, d2, ys, N2);
		else
			s_g_d(tc, tymin, tmax, tymax, jmin, jmax, t20, u2, diff2, b2, d2, ys, N2);
	}

	// (2) and (19)
	const double pt = ((std::min(t10, t20) + tc) / 2.);

	// (26)
	temp1 = voxel_index(pt, diff1, d1, apu_tx);
	// (27)
	temp2 = voxel_index(pt, diff2, d2, apu_ty);

	// (28)
	t1u = d1 / fabs(diff1);
	t2u = d2 / fabs(diff2);

	return false;
}

// Compute whether the ray intercepts the FOV (TYPE == 0), compute the voxel indices of the first voxel intercepted, 
// total number of voxels traversed and coefficients needed for Siddon's algorithm (3D)
bool siddon_pre_loop_3D(const double bx, const double by, const double bz, const double x_diff, const double y_diff, const double z_diff, 
	const double maxxx, const double maxyy, const double bzb, const double dx, const double dy, const double dz, 
	const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, int32_t& tempi, int32_t& tempj, int32_t& tempk, double& tyu, double& txu, double& tzu, 
	uint32_t& Np, const int TYPE, const Det detectors, double& tc, int32_t& iu, int32_t& ju, int32_t& ku, double& tx0, double& ty0, double& tz0) {

	const double apu_tx = bx - detectors.xs;
	const double apu_ty = by - detectors.ys;
	const double apu_tz = bz - detectors.zs;
	tx0 = (apu_tx) / (x_diff);
	ty0 = (apu_ty) / (y_diff);
	tz0 = (apu_tz) / (z_diff);
	const double txback = (maxxx - detectors.xs) / (x_diff);
	const double tyback = (maxyy - detectors.ys) / (y_diff);
	const double tzback = (bzb - detectors.zs) / (z_diff);

	const double txmin = std::min(tx0, txback);
	const double txmax = std::max(tx0, txback);
	const double tymin = std::min(ty0, tyback);
	const double tymax = std::max(ty0, tyback);
	const double tzmin = std::min(tz0, tzback);
	const double tzmax = std::max(tz0, tzback);

	tc = std::max(std::max(txmin, tzmin), tymin);
	const double tmax = std::min(std::min(txmax, tzmax), tymax);

	uint32_t imin, imax, jmin, jmax, kmin, kmax;

	if (TYPE == 0) {
		if (tc >= tmax) {
			return true;
		}

		if (detectors.xs < detectors.xd)
			d_g_s_precomp(tc, txmin, tmax, txmax, imin, imax, tx0, iu, x_diff, bx, dx, detectors.xs, Nx);
		else
			s_g_d_precomp(tc, txmin, tmax, txmax, imin, imax, tx0, iu, x_diff, bx, dx, detectors.xs, Nx);

		if (detectors.ys < detectors.yd)
			d_g_s_precomp(tc, tymin, tmax, tymax, jmin, jmax, ty0, ju, y_diff, by, dy, detectors.ys, Ny);
		else
			s_g_d_precomp(tc, tymin, tmax, tymax, jmin, jmax, ty0, ju, y_diff, by, dy, detectors.ys, Ny);

		if (detectors.zs < detectors.zd)
			d_g_s_precomp(tc, tzmin, tmax, tzmax, kmin, kmax, tz0, ku, z_diff, bz, dz, detectors.zs, Nz);
		else
			s_g_d_precomp(tc, tzmin, tmax, tzmax, kmin, kmax, tz0, ku, z_diff, bz, dz, detectors.zs, Nz);

		Np = (kmax - kmin + 1) + (jmax - jmin + 1) + (imax - imin + 1);
	}
	else {
		if (detectors.xs < detectors.xd)
			d_g_s(tc, txmin, tmax, txmax, imin, imax, tx0, iu, x_diff, bx, dx, detectors.xs, Nx);
		else
			s_g_d(tc, txmin, tmax, txmax, imin, imax, tx0, iu, x_diff, bx, dx, detectors.xs, Nx);

		if (detectors.ys < detectors.yd)
			d_g_s(tc, tymin, tmax, tymax, jmin, jmax, ty0, ju, y_diff, by, dy, detectors.ys, Ny);
		else
			s_g_d(tc, tymin, tmax, tymax, jmin, jmax, ty0, ju, y_diff, by, dy, detectors.ys, Ny);

		if (detectors.zs < detectors.zd)
			d_g_s(tc, tzmin, tmax, tzmax, kmin, kmax, tz0, ku, z_diff, bz, dz, detectors.zs, Nz);
		else
			s_g_d(tc, tzmin, tmax, tzmax, kmin, kmax, tz0, ku, z_diff, bz, dz, detectors.zs, Nz);
	}

	const double pt = ((std::min(std::min(tz0, ty0), tx0) + tc) / 2.);

	tempi = voxel_index(pt, x_diff, dx, apu_tx);
	tempj = voxel_index(pt, y_diff, dy, apu_ty);
	tempk = voxel_index(pt, z_diff, dz, apu_tz);

	txu = dx / fabs(x_diff);
	tyu = dy / fabs(y_diff);
	tzu = dz / fabs(z_diff);

	return false;
}

// Compute the total distance (and optionally forward projection) for the orthogonal ray (2D)
void orth_distance_full(const int32_t tempi, const uint32_t Nx, const double y_diff, const double x_diff, const double y_center, const double* x_center, 
	const double kerroin, const double length_, double& temp, const uint32_t tempijk, const uint32_t NN, const int32_t tempj, 
	const double local_sino, double& ax, const double* osem_apu, const bool no_norm, const bool RHS, const bool SUMMA, const bool OMP, const bool PRECOMP, 
	double* rhs, double* Summ, mwIndex* indices, std::vector<double>& elements, std::vector<uint32_t>& v_indices, size_t& idx, uint64_t N2) {

	const double diff = kerroin - x_diff * y_center;
	for (int32_t uu = tempi; uu >= 0; uu--) {
		double local_ele = compute_element_orth_mfree(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		const uint32_t local_ind = compute_ind_orth_mfree(static_cast<uint32_t>(uu), tempijk, NN);
		if (RHS) {
			local_ele *= temp;
#pragma omp atomic
			rhs[local_ind] += (local_ele * ax);
			if (no_norm == false) {
#pragma omp atomic
				Summ[local_ind] += local_ele;
			}
		}
		else if (SUMMA) {
			local_ele *= temp;
#pragma omp atomic
			Summ[local_ind] += local_ele;
		}
		else {
			temp += local_ele;
			if (OMP) {
				if (local_sino > 0.) {
					denominator_mfree(local_ele, ax, osem_apu[local_ind]);
				}
			}
			else if (PRECOMP) {
				rhs[N2 + idx] = (local_ele);
				indices[N2 + idx] = local_ind;
				idx++;
			}
			else {
				elements.emplace_back(local_ele);
				v_indices.emplace_back(local_ind);
				idx++;
			}
		}
	}
	for (uint32_t uu = static_cast<uint32_t>(tempi) + 1u; uu < Nx; uu++) {
		double local_ele = compute_element_orth_mfree(diff, y_diff, x_center[uu], length_);
		if (local_ele <= THR)
			break;
		const uint32_t local_ind = compute_ind_orth_mfree(uu, tempijk, NN);
		if (RHS) {
			local_ele *= temp;
#pragma omp atomic
			rhs[local_ind] += (local_ele * ax);
			if (no_norm == false) {
#pragma omp atomic
				Summ[local_ind] += local_ele;
			}
		}
		else if (SUMMA) {
			local_ele *= temp;
#pragma omp atomic
			Summ[local_ind] += local_ele;
		}
		else {
			temp += local_ele;
			if (OMP) {
				if (local_sino > 0.) {
					denominator_mfree(local_ele, ax, osem_apu[local_ind]);
				}
			}
			else if (PRECOMP) {
				rhs[N2 + idx] = (local_ele);
				indices[N2 + idx] = local_ind;
				idx++;
			}
			else {
				elements.emplace_back(local_ele);
				v_indices.emplace_back(local_ind);
				idx++;
			}
		}
	}
}

// Compute the total distance (and optionally forward projection) for the orthogonal ray (3D)
void orth_distance_3D_full(const int32_t tempi, const uint32_t Nx, const uint32_t Nz, const double y_diff, const double x_diff, const double z_diff,
	const double y_center, const double* x_center, const double* z_center, double& temp, const uint32_t tempijk, const uint32_t NN,
	const int32_t tempj, int32_t tempk, const double local_sino, double& ax, const double* osem_apu,
	const Det detectors, const uint32_t Nyx, const double crystal_size_z, const int32_t dec, const int32_t iu, 
	const bool no_norm, const bool RHS, const bool SUMMA, const bool OMP, const bool PRECOMP, double* rhs, double* Summ, mwIndex* indices, 
	std::vector<double>& elements, std::vector<uint32_t>& v_indices, size_t& idx, uint64_t N2) {

	bool loppu1 = true;
	bool loppu2 = true;
	uint8_t loppu3 = 0u;
	uint8_t loppu4 = 0u;
	bool pass1 = true;
	bool pass = false;
	int32_t alku = tempk;
	int32_t alku2 = tempk;
	int32_t alku3 = tempk;
	int32_t alku4 = tempk + 1;
	for (int32_t uu = tempi; uu >= 0; uu--) {
		//for (int32_t zz = tempk; zz >= 0; zz--) {
		for (int32_t zz = alku; zz >= 0; zz--) {
			double local_ele = compute_element_orth_3D(detectors, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
			if (local_ele <= THR) {
				if (loppu1 && zz >= alku3 - dec && !(pass && alku3 == tempk)) {
					continue;
				}
				else {
					if (loppu1 && loppu3 == 1u && loppu4 == 0u)
						pass1 = false;
					break;
				}
			}
			const uint32_t local_ind = compute_ind_orth_mfree_3D(static_cast<uint32_t>(uu), tempijk, static_cast<uint32_t>(zz), NN, Nyx);
			if (RHS) {
				local_ele *= temp;
#pragma omp atomic
				rhs[local_ind] += (local_ele * ax);
				if (no_norm == false) {
#pragma omp atomic
					Summ[local_ind] += local_ele;
				}
			}
			else if (SUMMA) {
				local_ele *= temp;
#pragma omp atomic
				Summ[local_ind] += local_ele;
			}
			else {
				temp += local_ele;
				if (OMP) {
					if (local_sino > 0.) {
						denominator_mfree(local_ele, ax, osem_apu[local_ind]);
					}
				}
				else if (PRECOMP) {
					rhs[N2 + idx] = (local_ele);
					indices[N2 + idx] = local_ind;
					idx++;
				}
				else {
					elements.emplace_back(local_ele);
					v_indices.emplace_back(local_ind);
					idx++;
				}
			}
			if (loppu1) {
				alku = zz;
				loppu1 = false;
				loppu3 = 1u;
				if (zz == tempk)
					pass = true;
			}
			alku3 = zz;
		}
		if (pass) {
			for (int32_t zz = alku2 + 1; zz < Nz; zz++) {
			//for (uint32_t zz = static_cast<uint32_t>(tempk) + 1u; zz < Nz; zz++) {
				double local_ele = compute_element_orth_3D(detectors, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
				if (local_ele <= THR) {
					if (loppu2 && zz <= alku4 + dec && (zz > (tempk + 1) || loppu1)) {
						continue;
					}
					else {
						if ((loppu2 && loppu4 == 1u) || zz == (tempk + 1)) {
							loppu4 = 0u;
							pass = false;
						}
						break;
					}
				}
				const uint32_t local_ind = compute_ind_orth_mfree_3D(static_cast<uint32_t>(uu), tempijk, static_cast<uint32_t>(zz), NN, Nyx);
				if (RHS) {
					local_ele *= temp;
#pragma omp atomic
					rhs[local_ind] += (local_ele * ax);
					if (no_norm == false) {
#pragma omp atomic
						Summ[local_ind] += local_ele;
					}
				}
				else if (SUMMA) {
					local_ele *= temp;
#pragma omp atomic
					Summ[local_ind] += local_ele;
				}
				else {
					temp += local_ele;
					if (OMP) {
						if (local_sino > 0.) {
							denominator_mfree(local_ele, ax, osem_apu[local_ind]);
						}
					}
					else if (PRECOMP) {
						rhs[N2 + idx] = (local_ele);
						indices[N2 + idx] = local_ind;
						idx++;
					}
					else {
						elements.emplace_back(local_ele);
						v_indices.emplace_back(local_ind);
						idx++;
					}
				}
				if (loppu2) {
					alku2 = zz - 1;
					loppu2 = false;
					loppu4 = 1u;
					pass1 = true;
				}
				alku4 = zz;
			}
		}
		loppu1 = true;
		loppu2 = true;
		if (pass1) {
			continue;
		}
		else {
			break;
		}
	}
	loppu3 = 0u;
	pass1 = true;
	pass = false;
	loppu4 = 0u;
	alku2 = tempk;
	alku = tempk;
	alku3 = tempk;
	alku4 = tempk + 1;
	for (uint32_t uu = static_cast<uint32_t>(tempi) + 1u; uu < Nx; uu++) {
		for (int32_t zz = alku; zz >= 0; zz--) {
		//for (int32_t zz = tempk; zz >= 0; zz--) {
			double local_ele = compute_element_orth_3D(detectors, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
			if (local_ele <= THR) {
				if (loppu1 && zz >= alku3 - dec && !(pass && alku3 == tempk)) {
					continue;
				}
				else {
					if (loppu1 && loppu3 == 1u && loppu4 == 0u)
						pass1 = false;
					break;
				}
			}
			const uint32_t local_ind = compute_ind_orth_mfree_3D(static_cast<uint32_t>(uu), tempijk, static_cast<uint32_t>(zz), NN, Nyx);
			if (RHS) {
				local_ele *= temp;
#pragma omp atomic
				rhs[local_ind] += (local_ele * ax);
				if (no_norm == false) {
#pragma omp atomic
					Summ[local_ind] += local_ele;
				}
			}
			else if (SUMMA) {
				local_ele *= temp;
#pragma omp atomic
				Summ[local_ind] += local_ele;
			}
			else {
				temp += local_ele;
				if (OMP) {
					if (local_sino > 0.) {
						denominator_mfree(local_ele, ax, osem_apu[local_ind]);
					}
				}
				else if (PRECOMP) {
					rhs[N2 + idx] = (local_ele);
					indices[N2 + idx] = local_ind;
					idx++;
				}
				else {
					elements.emplace_back(local_ele);
					v_indices.emplace_back(local_ind);
					idx++;
				}
			}
			if (loppu1) {
				alku = zz;
				loppu1 = false;
				loppu3 = 1u;
				if (zz == tempk)
					pass = true;
			}
			alku3 = zz;
		}
		if (pass) {
			for (int32_t zz = alku2 + 1; zz < Nz; zz++) {
			//for (uint32_t zz = static_cast<uint32_t>(tempk) + 1u; zz < Nz; zz++) {
				double local_ele = compute_element_orth_3D(detectors, x_diff, y_diff, z_diff, crystal_size_z, x_center[uu], y_center, z_center[zz]);
				if (local_ele <= THR) {
					if (loppu2 && zz <= alku4 + dec && (zz > (tempk + 1) || loppu1)) {
						continue;
					}
					else {
						if ((loppu2 && loppu4 == 1u) || zz == (tempk + 1)) {
							loppu4 = 0u;
							pass = false;
						}
						break;
					}
				}
				const uint32_t local_ind = compute_ind_orth_mfree_3D(static_cast<uint32_t>(uu), tempijk, static_cast<uint32_t>(zz), NN, Nyx);
				if (RHS) {
					local_ele *= temp;
#pragma omp atomic
					rhs[local_ind] += (local_ele * ax);
					if (no_norm == false) {
#pragma omp atomic
						Summ[local_ind] += local_ele;
					}
				}
				else if (SUMMA) {
					local_ele *= temp;
#pragma omp atomic
					Summ[local_ind] += local_ele;
				}
				else {
					temp += local_ele;
					if (OMP) {
						if (local_sino > 0.) {
							denominator_mfree(local_ele, ax, osem_apu[local_ind]);
						}
					}
					else if (PRECOMP) {
						rhs[N2 + idx] = (local_ele);
						indices[N2 + idx] = local_ind;
						idx++;
					}
					else {
						elements.emplace_back(local_ele);
						v_indices.emplace_back(local_ind);
						idx++;
					}
				}
				if (loppu2) {
					alku2 = zz - 1;
					loppu2 = false;
					loppu4 = 1u;
					pass1 = true;
				}
			}
		}
		loppu1 = true;
		loppu2 = true;
		if (pass1) {
			continue;
		}
		else {
			break;
		}
	}
}


// Compute the nominator (backprojection)
void nominator_mfree(double& ax, const double Sino, const double epps, const double temp, const bool randoms_correction, const double* randoms, 
	const uint32_t lo) {
	if (ax == 0.)
		ax = epps;
	else
		ax *= temp;
	if (randoms_correction)
		ax += randoms[lo];
	ax = Sino / ax;
}

// Compute the current system matrix element in the orthogonal ray tracer (2D)
double compute_element_orth_mfree(const double x_diff, const double y_diff, const double x_center, const double length_) {
	const double local_ele = 1. - fabs(x_diff + y_diff * x_center) / length_;
	return local_ele;
}

// Compute the current system matrix index in the orthogonal ray tracer (2D)
uint32_t compute_ind_orth_mfree(const uint32_t tempi, const uint32_t tempijk, const uint32_t d_N) {
	uint32_t local_ind = tempi * d_N + tempijk;
	return local_ind;
}

// Compute the current system matrix index in the orthogonal ray tracer (3D)
uint32_t compute_ind_orth_mfree_3D(const uint32_t tempi, const uint32_t tempijk, const int tempk, const uint32_t d_N, const uint32_t Nyx) {
	uint32_t local_ind = tempi * d_N + tempijk + tempk * Nyx;
	return local_ind;
}

// Denominator (forward projection) for the type 4
void denominator_mfree(const double local_ele, double& axOSEM, const double d_OSEM) {
	axOSEM += (local_ele * d_OSEM);
}

// Compute the starting index in a perpendicular detector case
uint32_t perpendicular_start(const double d_b, const double d, const double d_d, const uint32_t d_N) {
	uint32_t tempi = 0u;
	double start = d_b - d + d_d;
	for (uint32_t ii = 0u; ii < d_N; ii++) {
		if (start > 0.) {
			tempi = ii;
			break;
		}
		start += d_d;
	}
	return tempi;
}

// Calculate the denominator (forward projection) in the perpendicular case in orthogonal ray tracer (2D case)
void orth_distance_denominator_perpendicular_mfree(const double diff2, const double* center1, const double kerroin,
	const double length_, double& temp, const uint32_t d_attenuation_correction, double& ax, const double d_b, const double d, const double d_d1,
	const uint32_t d_N1, const uint32_t d_N2, const uint32_t z_loop, const double* d_atten, const double local_sino, const uint32_t d_N, const uint32_t d_NN,
	const double* d_OSEM) {

	const uint32_t zz = z_loop * d_N2 * d_N1;
	const uint32_t apu = perpendicular_start(d_b, d, d_d1, d_N1);
	double jelppi = 0.;
	for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
		const double local_ele = compute_element_orth_mfree(kerroin, diff2, center1[uu], length_);
		if (local_ele <= THR)
			break;
		temp += (local_ele * static_cast<double>(d_N2));
		uint32_t local_ind = uu * d_N + zz;
		for (uint32_t kk = 0u; kk < d_N2; kk++) {
			if (d_attenuation_correction && uu == static_cast<int32_t>(apu))
				jelppi += (d_d1 * -d_atten[local_ind]);
			if (local_sino > 0.f) {
				denominator_mfree(local_ele, ax, d_OSEM[local_ind]);
			}
			local_ind += d_NN;
		}
	}
	for (uint32_t hh = apu + 1; hh < d_N1; hh++) {
		const double local_ele = compute_element_orth_mfree(kerroin, diff2, center1[hh], length_);
		if (local_ele <= THR)
			break;
		temp += (local_ele * static_cast<double>(d_N2));
		uint32_t local_ind = hh * d_N + zz;
		if (local_sino > 0.f) {
			for (uint32_t kk = 0u; kk < d_N2; kk++) {
				denominator_mfree(local_ele, ax, d_OSEM[local_ind]);
				local_ind += d_NN;
			}
		}
	}
	temp = 1. / temp;
	if (d_attenuation_correction)
		temp *= jelppi;
}

// Calculate the normalization factor and RHS in the perpendicular case in orthogonal ray tracer (2D case)
void orth_distance_rhs_perpendicular_mfree(const double diff2, const double* center1, const double kerroin,
	const double length_, const double temp, double& ax, const double d_b, const double d, const double d_d1, const uint32_t d_N1, const uint32_t d_N2,
	const uint32_t z_loop, const uint32_t d_N, const uint32_t d_NN, const bool no_norm, double* rhs, double* Summ) {

	const uint32_t zz = z_loop * d_N2 * d_N1;
	const uint32_t apu = perpendicular_start(d_b, d, d_d1, d_N1);
	for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
		double local_ele = compute_element_orth_mfree(kerroin, diff2, center1[uu], length_);
		if (local_ele <= THR)
			break;
		uint32_t local_ind = uu * d_N + zz;
		local_ele *= temp;
		for (uint32_t kk = 0u; kk < d_N2; kk++) {
#pragma omp atomic
			rhs[local_ind] += (local_ele * ax);
			if (no_norm == 0) {
#pragma omp atomic
				Summ[local_ind] += local_ele;
			}
			local_ind += d_NN;
		}
	}
	for (uint32_t uu = apu + 1; uu < d_N1; uu++) {
		double local_ele = compute_element_orth_mfree(kerroin, diff2, center1[uu], length_);
		if (local_ele <= THR)
			break;
		uint32_t local_ind = uu * d_N + zz;
		local_ele *= temp;
		for (uint32_t kk = 0u; kk < d_N2; kk++) {
			if (no_norm == 0) {
#pragma omp atomic
				Summ[local_ind] += local_ele;
			}
#pragma omp atomic
			rhs[local_ind] += (local_ele * ax);
			local_ind += d_NN;
		}
	}
}

// Calculate the normalization factor in the perpendicular case in orthogonal ray tracer (2D case)
void orth_distance_summ_perpendicular_mfree(const double diff2, const double* center1, const double kerroin, const double length_, const double temp,
	double ax, const double d_b, const double d, const double d_d1, const uint32_t d_N1, const uint32_t d_N2, const uint32_t z_loop, const uint32_t d_N, 
	const uint32_t d_NN, double* Summ) {

	const uint32_t zz = z_loop * d_N2 * d_N1;
	const uint32_t apu = perpendicular_start(d_b, d, d_d1, d_N1);
	for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
		double local_ele = compute_element_orth_mfree(kerroin, diff2, center1[uu], length_);
		if (local_ele <= THR)
			break;
		uint32_t local_ind = uu * d_N + zz;
		local_ele *= temp;
		for (uint32_t kk = 0u; kk < d_N2; kk++) {
#pragma omp atomic
			Summ[local_ind] += local_ele;
			local_ind += d_NN;
		}
	}
	for (uint32_t uu = apu + 1; uu < d_N1; uu++) {
		double local_ele = compute_element_orth_mfree(kerroin, diff2, center1[uu], length_);
		if (local_ele <= THR)
			break;
		uint32_t local_ind = uu * d_N + zz;
		local_ele *= temp;
		for (uint32_t kk = 0u; kk < d_N2; kk++) {
#pragma omp atomic
			Summ[local_ind] += local_ele;
			local_ind += d_NN;
		}
	}
}

// Calculate the denominator (forward projection) in the perpendicular case in orthogonal ray tracer (3D case)
void orth_distance_denominator_perpendicular_mfree_3D(const double* center1, const double center2, const double* z_center, double& temp, 
	const uint32_t d_attenuation_correction, double& ax, const double d_b, const double d, const double d_d1, const uint32_t d_N1, 
	const uint32_t d_N2, const uint32_t z_loop, const double* d_atten, const double local_sino, const uint32_t d_N, const uint32_t d_NN, 
	const double* d_OSEM, Det detectors, const double xl, const double yl, const double zl, const double crystal_size_z, const uint32_t Nyx, 
	const uint32_t Nz) {

	//const uint32_t zz = z_loop * d_N2 * d_N1;
	const uint32_t apu = perpendicular_start(d_b, d, d_d1, d_N1);
	double jelppi = 0.;
	for (int32_t zz = static_cast<int32_t>(z_loop); zz >= 0; zz--) {
		for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
			double local_ele = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			temp += (local_ele * d_N2);
			uint32_t local_ind = uu * d_N + zz * Nyx;
			for (uint32_t kk = 0u; kk < d_N2; kk++) {
				if (d_attenuation_correction && uu == static_cast<int32_t>(apu) && zz == static_cast<int32_t>(z_loop))
					jelppi += (d_d1 * -d_atten[local_ind]);
				if (local_sino > 0.f) {
					denominator_mfree(local_ele, ax, d_OSEM[local_ind]);
				}
				local_ind += d_NN;
			}
		}
		for (uint32_t uu = apu + 1; uu < d_N1; uu++) {
			double local_ele = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			temp += (local_ele * d_N2);
			uint32_t local_ind = uu * d_N + zz * Nyx;
			if (local_sino > 0.f) {
				for (uint32_t kk = 0u; kk < d_N2; kk++) {
					denominator_mfree(local_ele, ax, d_OSEM[local_ind]);
					local_ind += d_NN;
				}
			}
		}
	}
	for (uint32_t zz = z_loop + 1u; zz < Nz; zz++) {
		for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
			double local_ele = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			temp += (local_ele * d_N2);
			uint32_t local_ind = uu * d_N + zz * Nyx;
			if (local_sino > 0.f) {
				for (uint32_t kk = 0u; kk < d_N2; kk++) {
					denominator_mfree(local_ele, ax, d_OSEM[local_ind]);
					local_ind += d_NN;
				}
			}
		}
		for (uint32_t uu = apu + 1; uu < d_N1; uu++) {
			double local_ele = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			temp += (local_ele * d_N2);
			uint32_t local_ind = uu * d_N + zz * Nyx;
			if (local_sino > 0.f) {
				for (uint32_t kk = 0u; kk < d_N2; kk++) {
					denominator_mfree(local_ele, ax, d_OSEM[local_ind]);
					local_ind += d_NN;
				}
			}
		}
	}
	temp = 1. / temp;
	if (d_attenuation_correction)
		temp *= jelppi;
}

// Calculate the normalization factor and RHS in the perpendicular case in orthogonal ray tracer (3D case)
void orth_distance_rhs_perpendicular_mfree_3D(const double* center1, const double center2, const double* z_center, const double temp, double& ax, 
	const double d_b, const double d, const double d_d1, const uint32_t d_N1, const uint32_t d_N2, const uint32_t z_loop, const uint32_t d_N, 
	const uint32_t d_NN, const bool no_norm, double* rhs, double* Summ, Det detectors, const double xl, const double yl, const double zl, 
	const double crystal_size_z, const uint32_t Nyx, const uint32_t Nz) {

	//const uint32_t zz = z_loop * d_N2 * d_N1;
	const uint32_t apu = perpendicular_start(d_b, d, d_d1, d_N1);
	for (int32_t zz = static_cast<int32_t>(z_loop); zz >= 0; zz--) {
		for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
			double local_ele = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint32_t local_ind = uu * d_N + zz * Nyx;
			local_ele *= temp;
			for (uint32_t kk = 0u; kk < d_N2; kk++) {
#pragma omp atomic
				rhs[local_ind] += (local_ele * ax);
				if (no_norm == 0) {
#pragma omp atomic
					Summ[local_ind] += local_ele;
				}
				local_ind += d_NN;
			}
		}
		for (uint32_t uu = apu + 1; uu < d_N1; uu++) {
			double local_ele = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint32_t local_ind = uu * d_N + zz * Nyx;
			local_ele *= temp;
			for (uint32_t kk = 0u; kk < d_N2; kk++) {
				if (no_norm == 0) {
#pragma omp atomic
					Summ[local_ind] += local_ele;
				}
#pragma omp atomic
				rhs[local_ind] += (local_ele * ax);
				local_ind += d_NN;
			}
		}
	}
	for (uint32_t zz = z_loop + 1u; zz < Nz; zz++) {
		for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
			double local_ele = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint32_t local_ind = uu * d_N + zz * Nyx;
			local_ele *= temp;
			for (uint32_t kk = 0u; kk < d_N2; kk++) {
#pragma omp atomic
				rhs[local_ind] += (local_ele * ax);
				if (no_norm == 0) {
#pragma omp atomic
					Summ[local_ind] += local_ele;
				}
				local_ind += d_NN;
			}
		}
		for (uint32_t uu = apu + 1; uu < d_N1; uu++) {
			double local_ele = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint32_t local_ind = uu * d_N + zz * Nyx;
			local_ele *= temp;
			for (uint32_t kk = 0u; kk < d_N2; kk++) {
				if (no_norm == 0) {
#pragma omp atomic
					Summ[local_ind] += local_ele;
				}
#pragma omp atomic
				rhs[local_ind] += (local_ele * ax);
				local_ind += d_NN;
			}
		}
	}
}

// Calculate the normalization factor in the perpendicular case in orthogonal ray tracer (3D case)
void orth_distance_summ_perpendicular_mfree_3D(const double* center1, const double center2, const double* z_center, const double temp,
	double ax, const double d_b, const double d, const double d_d1, const uint32_t d_N1, const uint32_t d_N2, const uint32_t z_loop, 
	const uint32_t d_N, const uint32_t d_NN, double* Summ, Det detectors, const double xl, const double yl, const double zl, const double crystal_size_z, 
	const uint32_t Nyx, const uint32_t Nz) {

	//const uint32_t zz = z_loop * d_N2 * d_N1;
	const uint32_t apu = perpendicular_start(d_b, d, d_d1, d_N1);
	for (int32_t zz = static_cast<int32_t>(z_loop); zz >= 0; zz--) {
		for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
			double local_ele = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint32_t local_ind = uu * d_N + zz * Nyx;
			local_ele *= temp;
			for (uint32_t kk = 0u; kk < d_N2; kk++) {
#pragma omp atomic
				Summ[local_ind] += local_ele;
				local_ind += d_NN;
			}
		}
		for (uint32_t uu = apu + 1; uu < d_N1; uu++) {
			double local_ele = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint32_t local_ind = uu * d_N + zz * Nyx;
			local_ele *= temp;
			for (uint32_t kk = 0u; kk < d_N2; kk++) {
#pragma omp atomic
				Summ[local_ind] += local_ele;
				local_ind += d_NN;
			}
		}
	}
	for (uint32_t zz = z_loop + 1u; zz < Nz; zz++) {
		for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
			double local_ele = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint32_t local_ind = uu * d_N + zz * Nyx;
			local_ele *= temp;
			for (uint32_t kk = 0u; kk < d_N2; kk++) {
#pragma omp atomic
				Summ[local_ind] += local_ele;
				local_ind += d_NN;
			}
		}
		for (uint32_t uu = apu + 1; uu < d_N1; uu++) {
			double local_ele = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
			uint32_t local_ind = uu * d_N + zz * Nyx;
			local_ele *= temp;
			for (uint32_t kk = 0u; kk < d_N2; kk++) {
#pragma omp atomic
				Summ[local_ind] += local_ele;
				local_ind += d_NN;
			}
		}
	}
}