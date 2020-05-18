/**************************************************************************
* Includes various functions required by the original Siddon, improved 
* Siddon, orthogonal distance-based and volume-based ray tracers. Brackets
* signify equation numbers from the reference below.
* 
* References: 
* [1] Jacobs, F., Sundermann, E., De Sutter, B., Christiaens, M. Lemahieu, 
* I. (1998). A Fast Algorithm to Calculate the Exact Radiological Path 
* through a Pixel or Voxel Space. Journal of computing and information 
* technology, 6 (1), 89-94.
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

// Compute the orthogonal distance between a point and a line (in 3D)
double compute_element_orth_3D(Det detectors, const double xl, const double yl, const double zl, const double crystal_size_z, 
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

	// Normalize the distance
	return (1. - norm(x1, y1, z1) / crystal_size_z);
}

// Gaussian weight
// Compute the orthogonal distance between a point and a line (in 3D)
//double compute_element_orth_3D(Det detectors, const double xl, const double yl, const double zl, const double crystal_size_z,
//	const double xp, const double yp, const double zp) {
//
//	double x1, y1, z1, x0, y0, z0;
//
//	x0 = xp - detectors.xs;
//	y0 = yp - detectors.ys;
//	z0 = zp - detectors.zs;
//
//	// Cross product
//	x1 = yl * z0 - zl * y0;
//	y1 = zl * x0 - xl * z0;
//	z1 = xl * y0 - yl * x0;
//
//	const double normi = norm(x1, y1, z1);
//
//	double gauss = 0.;
//
//	if (normi < crystal_size_z)
//		gauss = (1. - std::exp(-(normi * normi) / (2. * (crystal_size_z / 2.35482) * (crystal_size_z / 2.35482))));
//
//	return gauss;
//}

// Compute the Euclidean norm of a vector
double norm(const double x, const double y, const double z) {
	return sqrt(x * x + y * y + z * z);
}

// Compute the current system matrix element in the orthogonal ray tracer (2D)
double compute_element_orth_mfree(const double x_diff, const double y_diff, const double x_center, const double length_) {
	const double local_ele = 1. - fabs(x_diff + y_diff * x_center) / length_;

	// Gaussian weight
	//const double normi = fabs(x_diff + y_diff * x_center);
	//const double local_ele = (1. - std::exp(-(normi * normi) / (2. * (length_ / 2.35482) * (length_ / 2.35482))));

	return local_ele;
}

// Compute the current system matrix index in the orthogonal ray tracer (2D)
uint32_t compute_ind_orth_mfree(const uint32_t tempi, const uint32_t tempijk, const uint32_t d_N) {
	uint32_t local_ind = tempi * d_N + tempijk;
	return local_ind;
}

// Compute the current system matrix index in the orthogonal ray tracer (3D)
uint32_t compute_ind_orth_mfree_3D(const uint32_t tempi, const uint32_t tempijk, const uint32_t tempk, const uint32_t d_N, const uint32_t Nyx) {
	uint32_t local_ind = tempi * d_N + tempijk + tempk * Nyx;
	return local_ind;
}

void computeIndices(const bool RHS, const bool SUMMA, const bool OMP, const bool PRECOMP, const bool DISCARD, double local_ele, double& temp, double& ax, 
	const bool no_norm, double* Summ, double* rhs, const double local_sino, const double* osem_apu, const uint64_t N2, mwIndex* indices, 
	std::vector<double>& elements, std::vector<uint32_t>& v_indices, size_t& idx, const uint32_t local_ind, const uint64_t N22) {
	// Compute the total probability for both backprojection and sensitivity image
	if (RHS) {
		local_ele *= temp;
#pragma omp atomic
		rhs[local_ind] += (local_ele * ax);
		if (no_norm == false) {
#pragma omp atomic
			Summ[local_ind] += local_ele;
		}
	}
	// Compute the total probability for sensitivity image only (no measurements)
	else if (SUMMA) {
		local_ele *= temp;
#pragma omp atomic
		Summ[local_ind] += local_ele;
	}
	else {
		// Preliminary pass to compute the total length of the ray in FOV and forward projection
		// Implementation 4
		if (OMP) {
			if (local_sino > 0.) {
				denominator_mfree(local_ele, ax, osem_apu[local_ind]);
			}
		}
		// Implementation 1
		else if (PRECOMP) {
			rhs[N2 + idx] = local_ele;
			indices[N2 + idx] = local_ind;
			idx++;
		}
		// Precompute phase
		else if (DISCARD)
			idx++;
		// Implementation 1 without precomputation
		else {
			elements.emplace_back(local_ele);
			v_indices.emplace_back(local_ind);
			idx++;
		}
		temp += local_ele;
	}
}

// Get the detector coordinates for the current raw list-mode measurement
void get_detector_coordinates_raw(const uint32_t det_per_ring, const double* x, const double* y, const double* z, Det& detectors,
	const uint16_t* L, const size_t ll, const uint32_t* pseudos, const uint32_t pRows, const bool list_mode_format) {
	uint32_t ps;
	if (list_mode_format) {
		detectors.zs = z[ll];
		detectors.zd = z[ll + det_per_ring];
		detectors.xs = x[ll];
		detectors.xd = x[ll + det_per_ring];
		detectors.ys = y[ll];
		detectors.yd = y[ll + det_per_ring];
	}
	else {
		const uint32_t detektorit1 = static_cast<uint32_t>(L[ll * 2ULL]) - 1u;
		const uint32_t detektorit2 = static_cast<uint32_t>(L[ll * 2ULL + 1ULL]) - 1u;

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

		detectors.xs = x[detektorit1 - det_per_ring * (loop1)];
		detectors.xd = x[detektorit2 - det_per_ring * (loop2)];
		detectors.ys = y[detektorit1 - det_per_ring * (loop1)];
		detectors.yd = y[detektorit2 - det_per_ring * (loop2)];
	}
}

// Get the detector coordinates for the current raw list-mode measurement (multi-ray)
void get_detector_coordinates_raw_N(const uint32_t det_per_ring, const double* x, const double* y, const double* z, Det& detectors,
	const uint16_t* L, const size_t ll, const uint32_t* pseudos, const uint32_t pRows, const uint16_t lor, const double cr_pz, 
	const uint16_t n_rays, const uint16_t n_rays3D) {
	uint32_t ps;
	const uint32_t detektorit1 = static_cast<uint32_t>(L[ll * 2ULL]) - 1u;
	const uint32_t detektorit2 = static_cast<uint32_t>(L[ll * 2ULL + 1ULL]) - 1u;

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

	if (n_rays3D > 1) {
		uint16_t rays3D = n_rays3D;
		if (n_rays3D % 2 == 0)
			rays3D++;
		std::vector<int> r(rays3D, 0.);
		std::iota(r.begin(), r.end() + 1, -rays3D / 2);
		if (n_rays3D % 2 == 0)
			r.erase(r.begin() + rays3D / 2);
		const double rr = static_cast<double>(r[(lor - 1) % n_rays3D]);
		detectors.zs += (cr_pz * rr);
		detectors.zd += (cr_pz * rr);
	}
	if (n_rays > 1) {
		const uint16_t ll = (lor - 1) / n_rays3D;
		detectors.xs = x[detektorit1 - det_per_ring * loop1 + det_per_ring * ll];
		detectors.xd = x[detektorit2 - det_per_ring * loop2 + det_per_ring * ll];
		detectors.ys = y[detektorit1 - det_per_ring * loop1 + det_per_ring * ll];
		detectors.yd = y[detektorit2 - det_per_ring * loop2 + det_per_ring * ll];
	}
	else {
		detectors.xs = x[detektorit1 - det_per_ring * (loop1)];
		detectors.xd = x[detektorit2 - det_per_ring * (loop2)];
		detectors.ys = y[detektorit1 - det_per_ring * (loop1)];
		detectors.yd = y[detektorit2 - det_per_ring * (loop2)];
	}

}

// Get the detector coordinates for the current sinogram bin, no precomputations performed beforehand
void get_detector_coordinates_noalloc(const double* x, const double* y, const double* z, const uint32_t size_x, Det& detectors, int& ll, 
	const uint32_t* index, int& lz, const uint32_t TotSinos, size_t oo) {
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
	const uint16_t* z_index, const uint32_t TotSinos, const size_t oo) {
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
	const uint16_t* z_index, const uint32_t TotSinos, const size_t oo, const uint16_t lor, const double cr_pz, const uint16_t n_rays, const uint16_t n_rays3D) {
	// Sinogram data
	// Multiple axial rays
	if (n_rays3D > 1) {
		uint16_t rays3D = n_rays3D;
		// Increment by one if even
		if (n_rays3D % 2 == 0)
			rays3D++;
		// Create a vector of size n_rays3D
		std::vector<int> r(rays3D, 0.);
		// Populate the vector with values from -x...x
		std::iota(r.begin(), r.end() + 1, -rays3D / 2);
		// If even, remove the center 0
		if (n_rays3D % 2 == 0)
			r.erase(r.begin() + rays3D / 2);
		// Select the r-value for the current ray number
		const double rr = static_cast<double>(r[(lor - 1) % n_rays3D]);
		// Add the previously computed distance between rays
		// If, e.g. 2 axial rays, then rr will be -1 and 1
		// With 3 rays -1, 0, 1
		// z-detector values are always in the center of the crystal
		detectors.zs = z[z_index[oo]] + cr_pz * rr;
		if (z_index[oo] >= TotSinos) {
			detectors.zd = z[z_index[oo] - TotSinos] + cr_pz * rr;
		}
		else {
			detectors.zd = z[z_index[oo] + TotSinos] + cr_pz * rr;
		}
	}
	// No multiple axial rays
	else {
		detectors.zs = z[z_index[oo]];
		if (z_index[oo] >= TotSinos) {
			detectors.zd = z[z_index[oo] - TotSinos];
		}
		else {
			detectors.zd = z[z_index[oo] + TotSinos];
		}
	}
	// Multiple transaxial rays
	if (n_rays > 1) {
		// Transaxial coordinates are precomputed for the multiray case
		// Based on current ray number, determine the correct location of the coordinates
		// Input coordinates are a 3D matrix where first two dimensions are the sinogram coordinates
		// and the third dimension is the number of transaxial rays
		const uint16_t ll = (lor - 1) / n_rays3D * 2;
		detectors.xs = x[xy_index[oo] + size_x * ll];
		detectors.ys = y[xy_index[oo] + size_x * ll];
		if (xy_index[oo] >= size_x) {
			detectors.xd = x[xy_index[oo] + size_x * (ll - 1)];
			detectors.yd = y[xy_index[oo] + size_x * (ll - 1)];
		}
		else {
			detectors.xd = x[xy_index[oo] + size_x * (ll + 1)];
			detectors.yd = y[xy_index[oo] + size_x * (ll + 1)];
		}
	}
	// No transaxial rays
	else {
		detectors.ys = y[xy_index[oo]];
		detectors.xs = x[xy_index[oo]];
		if (xy_index[oo] >= size_x) {
			detectors.xd = x[xy_index[oo] - size_x];
			detectors.yd = y[xy_index[oo] - size_x];
		}
		else {
			detectors.xd = x[xy_index[oo] + size_x];
			detectors.yd = y[xy_index[oo] + size_x];
		}
	}
}

// Get the current ring in axial direction
uint32_t z_ring(const double zmax, const double zs, const double NSlices) {
	return static_cast<uint32_t>((zs / zmax) * (NSlices - 1.));
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

// Compute the sum in the attenuation correction (distance times attenuation coefficient)
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
	uint32_t& tempk, const uint32_t NN, const size_t lo, const double global_factor) {
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
	temp *= global_factor;


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
	// Voxel number
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
void orth_perpendicular_precompute(const uint32_t N1, const uint32_t N2, const double dd, const std::vector<double> vec, 
	const double* center1, const double center2, const double* z_center, const double kerroin, size_t &temp_koko, const Det detectors,
	const double xl, const double yl, const double zl, const int32_t tempk) {
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
		const double d_ort = compute_element_orth_3D(detectors, xl, yl, zl, kerroin, center1[uu], center2, z_center[tempk]);
		if (d_ort <= THR)
			break;
		koko1++;
	}
	temp_koko += koko1 * N1;
	for (uint32_t uu = apu + 1u; uu < N2; uu++) {
		const double d_ort = compute_element_orth_3D(detectors, xl, yl, zl, kerroin, center1[uu], center2, z_center[tempk]);
		if (d_ort <= THR)
			break;
		koko2++;
	}
	temp_koko += koko2 * N1;
}

// compute the orthogonal distance, perpendicular detectors
void orth_perpendicular_precompute_3D(const uint32_t N1, const uint32_t N2, const uint32_t Nz, const double dd, const std::vector<double> vec, 
	const double* center1, const double center2, const double* z_center, const double crystal_size_z, size_t& temp_koko, const Det detectors,
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
	for (int32_t uu = tempi + 1; uu < Nx; uu++) {
		const double d_ort = compute_element_orth_mfree(diff, y_diff, x_center[uu], length_);
		if (d_ort <= THR)
			break;
		temp_koko++;
	}
}

void orth_distance_precompute_3D(const int32_t tempi, const uint32_t N, const uint32_t Nz, const double y_diff, const double x_diff, const double z_diff, 
	const double* x_center, const double y_center, const double* z_center, const double crystal_size_z, const Det detectors, 
	uint16_t& temp_koko, const int32_t tempk, const size_t lo, const int32_t n_tempk, const int32_t dec, const int32_t decx) {

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
	for (int32_t uu = tempi + 1; uu < N; uu++) {
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
	//if (temp2 > N2) {
	//	mexPrintf("pt = %f\n", pt);
	//	mexPrintf("diff2 = %f\n", diff2);
	//	mexPrintf("d2 = %f\n", d2);
	//	mexPrintf("apu_ty = %f\n", apu_ty);
	//}

	if (TYPE == 0) {
		if (temp1 < 0 || temp1 >= N1 || temp2 < 0 || temp2 >= N2)
			return true;
	}

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

	if (TYPE == 0) {
		if (tempi < 0 || tempi >= Nx || tempj < 0 || tempj >= Ny || tempk < 0 || tempk >= Nz)
			return true;
	}

	txu = dx / fabs(x_diff);
	tyu = dy / fabs(y_diff);
	tzu = dz / fabs(z_diff);

	return false;
}


// Compute the total distance (and optionally forward projection) for the orthogonal ray (3D)
void orth_distance_3D_full(int32_t tempi, const uint32_t Nx, const uint32_t Nz, const double y_diff, const double x_diff, const double z_diff,
	const double* y_center, const double* x_center, const double* z_center, double& temp, const uint32_t NN, int32_t tempj, int32_t tempk, 
	const double local_sino, double& ax, const double* osem_apu, const Det detectors, const uint32_t Nyx, const double kerroin, 
	const bool no_norm, const bool RHS, const bool SUMMA, const bool OMP, const bool PRECOMP, const bool DISCARD, double* rhs, double* Summ, mwIndex* indices,
	std::vector<double>& elements, std::vector<uint32_t>& v_indices, size_t& idx, const uint32_t Ny, const uint32_t N1, const int start,
	const int32_t iu, const int32_t ju, const int loppu, std::vector<double>& store_elements, std::vector<uint32_t>& store_indices,
	const uint32_t tid, uint32_t& ind, uint64_t N2, uint64_t N22) {

	// Use the stored intersection lengths and voxel indices to compute the final probabilities
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
		// Increasing axial case
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
				double prev_local = 1.;
				for (xx = alku_x1; xx < Nx; xx++) {
					// Compute the normalized orthogonal distance <= 1
					double local_ele = compute_element_orth_3D(detectors, x_diff, y_diff, z_diff, kerroin, x_center[xx], y_center[yy1], z_center[zz]);
					// Always compute at least one additional orthogonal distance
					// determine if the voxel is further away from the ray than the previous one
					// Break the loop if the voxel distance increases
					if (local_ele <= THR && incr > 0 && prev_local > local_ele) {
						if (xx == alku_x1 + 1) {
							breikki1 = true;
						}
						break;
						//continue;
					}
					// This is computed at least once
					// Otherwise if voxels are getting closer to the ray
					else if (local_ele <= THR) {
						incr++;
						prev_local = local_ele;
						continue;
					}
					incr = 1;
					prev_local = local_ele;
					// Compute the linear index of the current voxel
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
				prev_local = 1.;
				for (xx = alku_x2; xx >= 0; xx--) {
					double local_ele = compute_element_orth_3D(detectors, x_diff, y_diff, z_diff, kerroin, x_center[xx], y_center[yy1], z_center[zz]);
					if (local_ele <= THR && incr > 0 && prev_local > local_ele) {
						if (xx == alku_x2 - 1) {
							breikki2 = true;
						}
						break;
						//continue;
					}
					else if (local_ele <= THR) {
						incr++;
						prev_local = local_ele;
						continue;
					}
					incr = 1;
					prev_local = local_ele;
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
						alku_x2 = alku_x1 - 1;
					}
					else {
						alku_x2 = uu2 + 1;
						alku_x1 = alku_x2 + 1;
					}
				}
				else {
					if (ju > 0) {
						alku_x2 = uu2 + 1;
						alku_x1 = alku_x2 + 1;
					}
					else {
						alku_x1 = uu1 - 1;
						alku_x2 = alku_x1 - 1;
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
				double prev_local = 1.;
				for (xx = alku_x1; xx < Nx; xx++) {
					double local_ele = compute_element_orth_3D(detectors, x_diff, y_diff, z_diff, kerroin, x_center[xx], y_center[yy2], z_center[zz]);
					if (local_ele <= THR && incr > 0 && prev_local > local_ele) {
						if (xx == alku_x1 + 1) {
							breikki1 = true;
						}
						break;
						//continue;
					}
					else if (local_ele <= THR) {
						incr++;
						prev_local = local_ele;
						continue;
					}
					incr = 1;
					prev_local = local_ele;
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
				prev_local = 1.;
				for (xx = alku_x2; xx >= 0; xx--) {
					double local_ele = compute_element_orth_3D(detectors, x_diff, y_diff, z_diff, kerroin, x_center[xx], y_center[yy2], z_center[zz]);
					if (local_ele <= THR && incr > 0 && prev_local > local_ele) {
						if (xx == alku_x2 - 1) {
							breikki2 = true;
						}
						break;
						//continue;
					}
					else if (local_ele <= THR) {
						incr++;
						prev_local = local_ele;
						continue;
					}
					incr = 1;
					prev_local = local_ele;
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
		//for (int zz = tempk - 1; zz >= loppu; zz--) {
		//	yy1 = 0;
		//	uu1 = 0;
		//	yy2 = 0;
		//	uu2 = 0;
		//	alku_y1 = tempj;
		//	alku_x1 = tempi;
		//	alku_y2 = tempj - 1;
		//	alku_x2 = tempi - 1;
		//	breikki1 = false;
		//	breikki2 = false;
		//	breikki3 = false;
		//	breikki4 = false;
		//	for (yy1 = alku_y1; yy1 < Ny; yy1++) {
		//		int xx = 0;
		//		int incr = 0;
		//		double prev_local = 1.;
		//		for (xx = alku_x1; xx < Nx; xx++) {
		//			double local_ele = compute_element_orth_3D(detectors, x_diff, y_diff, z_diff, kerroin, x_center[xx], y_center[yy1], z_center[zz]);
		//			if (local_ele <= THR && incr > 0 && prev_local > local_ele) {
		//				if (xx == alku_x1 + 1) {
		//					breikki1 = true;
		//				}
		//				//break;
		//				continue;
		//			}
		//			else if (local_ele <= THR) {
		//				incr++;
		//				prev_local = local_ele;
		//				continue;
		//			}
		//			incr = 1;
		//			prev_local = local_ele;
		//			const uint32_t local_ind = compute_ind_orth_mfree_3D(static_cast<uint32_t>(xx), yy1 * N1, static_cast<uint32_t>(zz), NN, Nyx);
		//			computeIndices(RHS, SUMMA, OMP, PRECOMP, DISCARD, local_ele, temp, ax, no_norm, Summ, rhs,
		//				local_sino, osem_apu, N2, indices, elements, v_indices, idx, local_ind, N22);
		//			if (!DISCARD && !PRECOMP) {
		//				store_elements[tid + ind] = local_ele;
		//				store_indices[tid + ind] = local_ind;
		//				ind++;
		//			}
		//		}
		//		uu1 = xx;
		//		xx = 0;
		//		incr = 0;
		//		prev_local = 1.;
		//		for (xx = alku_x2; xx >= 0; xx--) {
		//			double local_ele = compute_element_orth_3D(detectors, x_diff, y_diff, z_diff, kerroin, x_center[xx], y_center[yy1], z_center[zz]);
		//			if (local_ele <= THR && incr > 0 && prev_local > local_ele) {
		//				if (xx == alku_x2 - 1) {
		//					breikki2 = true;
		//				}
		//				//break;
		//				continue;
		//			}
		//			else if (local_ele <= THR) {
		//				incr++;
		//				prev_local = local_ele;
		//				continue;
		//			}
		//			incr = 1;
		//			prev_local = local_ele;
		//			const uint32_t local_ind = compute_ind_orth_mfree_3D(static_cast<uint32_t>(xx), yy1 * N1, static_cast<uint32_t>(zz), NN, Nyx);
		//			computeIndices(RHS, SUMMA, OMP, PRECOMP, DISCARD, local_ele, temp, ax, no_norm, Summ, rhs,
		//				local_sino, osem_apu, N2, indices, elements, v_indices, idx, local_ind, N22);
		//			if (!DISCARD && !PRECOMP) {
		//				store_elements[tid + ind] = local_ele;
		//				store_indices[tid + ind] = local_ind;
		//				ind++;
		//			}
		//		}
		//		uu2 = xx;
		//		if (iu > 0) {
		//			if (ju > 0) {
		//				alku_x1 = uu1 - 1;
		//				alku_x2 = uu1 - 2;
		//			}
		//			else {
		//				alku_x2 = uu2 + 1;
		//				alku_x1 = uu2 + 2;
		//			}
		//		}
		//		else {
		//			if (ju > 0) {
		//				alku_x2 = uu2 + 1;
		//				alku_x1 = uu2 + 2;
		//			}
		//			else {
		//				alku_x1 = uu1 - 1;
		//				alku_x2 = uu1 - 2;
		//			}
		//		}
		//		if (breikki1 && breikki2) {
		//			if (yy1 == alku_y1)
		//				breikki3 = true;
		//			break;
		//		}
		//	}
		//	breikki1 = false;
		//	breikki2 = false;
		//	alku_x1 = tempi;
		//	alku_x2 = tempi - 1;
		//	for (yy2 = alku_y2; yy2 >= 0; yy2--) {
		//		int xx = 0;
		//		int incr = 0;
		//		double prev_local = 1.;
		//		for (xx = alku_x1; xx < Nx; xx++) {
		//			double local_ele = compute_element_orth_3D(detectors, x_diff, y_diff, z_diff, kerroin, x_center[xx], y_center[yy2], z_center[zz]);
		//			if (local_ele <= THR && incr > 0 && prev_local > local_ele) {
		//				if (xx == alku_x1 + 1) {
		//					breikki1 = true;
		//				}
		//				//break;
		//				continue;
		//			}
		//			else if (local_ele <= THR) {
		//				incr++;
		//				prev_local = local_ele;
		//				continue;
		//			}
		//			incr = 1;
		//			prev_local = local_ele;
		//			const uint32_t local_ind = compute_ind_orth_mfree_3D(static_cast<uint32_t>(xx), yy2 * N1, static_cast<uint32_t>(zz), NN, Nyx);
		//			computeIndices(RHS, SUMMA, OMP, PRECOMP, DISCARD, local_ele, temp, ax, no_norm, Summ, rhs,
		//				local_sino, osem_apu, N2, indices, elements, v_indices, idx, local_ind, N22);
		//			if (!DISCARD && !PRECOMP) {
		//				store_elements[tid + ind] = local_ele;
		//				store_indices[tid + ind] = local_ind;
		//				ind++;
		//			}
		//		}
		//		uu1 = xx;
		//		xx = 0;
		//		incr = 0;
		//		prev_local = 1.;
		//		for (xx = alku_x2; xx >= 0; xx--) {
		//			double local_ele = compute_element_orth_3D(detectors, x_diff, y_diff, z_diff, kerroin, x_center[xx], y_center[yy2], z_center[zz]);
		//			if (local_ele <= THR && incr > 0 && prev_local > local_ele) {
		//				if (xx == alku_x2 - 1) {
		//					breikki2 = true;
		//				}
		//				//break;
		//				continue;
		//			}
		//			else if (local_ele <= THR) {
		//				incr++;
		//				prev_local = local_ele;
		//				continue;
		//			}
		//			incr = 1;
		//			prev_local = local_ele;
		//			const uint32_t local_ind = compute_ind_orth_mfree_3D(static_cast<uint32_t>(xx), yy2 * N1, static_cast<uint32_t>(zz), NN, Nyx);
		//			computeIndices(RHS, SUMMA, OMP, PRECOMP, DISCARD, local_ele, temp, ax, no_norm, Summ, rhs,
		//				local_sino, osem_apu, N2, indices, elements, v_indices, idx, local_ind, N22);
		//			if (!DISCARD && !PRECOMP) {
		//				store_elements[tid + ind] = local_ele;
		//				store_indices[tid + ind] = local_ind;
		//				ind++;
		//			}
		//		}
		//		uu2 = xx;
		//		if (iu > 0) {
		//			if (ju > 0) {
		//				alku_x2 = uu2 + 1;
		//				alku_x1 = uu2 + 2;
		//			}
		//			else {
		//				alku_x1 = uu1 - 1;
		//				alku_x2 = uu1 - 2;
		//			}
		//		}
		//		else {
		//			if (ju > 0) {
		//				alku_x1 = uu1 - 1;
		//				alku_x2 = uu1 - 2;
		//			}
		//			else {
		//				alku_x2 = uu2 + 1;
		//				alku_x1 = uu2 + 2;
		//			}
		//		}
		//		if (breikki1 && breikki2) {
		//			if (yy2 == alku_y2)
		//				breikki4 = true;
		//			break;
		//		}
		//	}
		//	if (breikki3 && breikki4) {
		//		break;
		//	}
		//}

	}
}


// Compute the nominator (backprojection)
void nominator_mfree(double& ax, const double Sino, const double epps, const double temp, const bool randoms_correction, const double* randoms, 
	const size_t lo) {
	if (ax == 0.)
		ax = epps;
	else
		ax *= temp;
	if (randoms_correction)
		ax += randoms[lo];
	ax = Sino / ax;
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
void orth_distance_denominator_perpendicular_mfree(const double* center1, const double center2, const double* z_center, const double kerroin,
	double& temp, const bool d_attenuation_correction, const bool normalization, double& ax, const double d_b, const double d, const double d_d1,
	const uint32_t d_N1, const uint32_t d_N2, const uint32_t z_loop, const double* d_atten, const double* norm_coef, const double local_sino, const uint32_t d_N, const uint32_t d_NN,
	const double* d_OSEM, const Det detectors, const double xl, const double yl, const double zl, std::vector<double>& store_elements, std::vector<uint32_t>& store_indices, 
	const uint32_t tid, uint32_t& ind, double* elements, mwIndex* indices, const size_t lo, const bool PRECOMPUTE, const double global_factor, const uint64_t N2) {

	const uint32_t zz = z_loop * d_N2 * d_N1;
	const uint32_t apu = perpendicular_start(d_b, d, d_d1, d_N1);
	double jelppi = 0.;
	for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
		double local_ele = compute_element_orth_3D(detectors, xl, yl, zl, kerroin, center1[uu], center2, z_center[z_loop]);
		if (local_ele <= THR)
			break;
		temp += (local_ele * static_cast<double>(d_N2));
		uint32_t local_ind = uu * d_N + zz;
		if (!PRECOMPUTE) {
			store_indices[tid + ind] = local_ind;
			store_elements[tid + ind] = local_ele;
			ind++;
		}
		for (uint32_t kk = 0u; kk < d_N2; kk++) {
			if (d_attenuation_correction && uu == static_cast<int32_t>(apu))
				jelppi += (d_d1 * -d_atten[local_ind]);
			if (local_sino > 0.) {
				denominator_mfree(local_ele, ax, d_OSEM[local_ind]);
			}
			local_ind += d_NN;
		}
	}
	for (uint32_t uu = apu + 1; uu < d_N1; uu++) {
		double local_ele = compute_element_orth_3D(detectors, xl, yl, zl, kerroin, center1[uu], center2, z_center[z_loop]);
		if (local_ele <= THR)
			break;
		temp += (local_ele * static_cast<double>(d_N2));
		uint32_t local_ind = uu * d_N + zz;
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
	temp = 1. / temp;
	if (d_attenuation_correction)
		temp *= exp(jelppi);
	if (normalization)
		temp *= norm_coef[lo];
	temp *= global_factor;
	if (PRECOMPUTE) {
		uint32_t hpk = N2;
		for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
			double local_ele = compute_element_orth_3D(detectors, xl, yl, zl, kerroin, center1[uu], center2, z_center[z_loop]);
			if (local_ele <= THR)
				break;
			local_ele *= temp;
			uint32_t local_ind = uu * d_N + zz;
			for (uint32_t kk = 0u; kk < d_N2; kk++) {
				indices[hpk] = local_ind;
				elements[hpk] = local_ele;
				local_ind += d_NN;
				hpk++;
			}
		}
		for (uint32_t uu = apu + 1; uu < d_N1; uu++) {
			double local_ele = compute_element_orth_3D(detectors, xl, yl, zl, kerroin, center1[uu], center2, z_center[z_loop]);
			if (local_ele <= THR)
				break;
			local_ele *= temp;
			uint32_t local_ind = uu * d_N + zz;
			for (uint32_t kk = 0u; kk < d_N2; kk++) {
				indices[hpk] = local_ind;
				elements[hpk] = local_ele;
				local_ind += d_NN;
				hpk++;
			}
		}
	}
}

// Calculate the normalization factor and RHS in the perpendicular case in orthogonal ray tracer (2D case)
void orth_distance_rhs_perpendicular_mfree(const double* center1, const double center2, const double* z_center, const double kerroin,
	const double temp, double& ax, const double d_b, const double d, const double d_d1, const uint32_t d_N1, const uint32_t d_N2,
	const uint32_t z_loop, const uint32_t d_N, const uint32_t d_NN, const bool no_norm, double* rhs, double* Summ, const bool RHS, const bool SUMMA, 
	const Det detectors, const double xl, const double yl, const double zl, const std::vector<double> store_elements, const std::vector<uint32_t> store_indices,
	const uint32_t tid, uint32_t ind, double* elements, mwIndex* indices, uint64_t N2) {

	const uint32_t zz = z_loop * d_N2 * d_N1;
	const uint32_t apu = perpendicular_start(d_b, d, d_d1, d_N1);
	uint32_t hpk = N2;
	//for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
	for (int32_t uu = 0; uu < ind; uu++) {
		double local_ele = store_elements[tid + uu];
		uint32_t local_ind = store_indices[tid + uu];
		local_ele *= temp;
		if (RHS) {
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
		else if (SUMMA) {
			for (uint32_t kk = 0u; kk < d_N2; kk++) {
#pragma omp atomic
				Summ[local_ind] += local_ele;
				local_ind += d_NN;
			}
		}
		else {
			for (uint32_t kk = 0u; kk < d_N2; kk++) {
				indices[hpk] = local_ind;
				elements[hpk] = local_ele;
				local_ind += d_NN;
				hpk++;
			}
		}
	}
}

// Calculate the denominator (forward projection) in the perpendicular case in orthogonal ray tracer (3D case)
void orth_distance_denominator_perpendicular_mfree_3D(const double* center1, const double center2, const double* z_center, double& temp, 
	const bool d_attenuation_correction, const bool normalization, double& ax, const double d_b, const double d, const double d_d1, const uint32_t d_N1,
	const uint32_t d_N2, const uint32_t z_loop, const double* d_atten, const double* norm_coef, const double local_sino, const uint32_t d_N, const uint32_t d_NN, 
	const double* d_OSEM, Det detectors, const double xl, const double yl, const double zl, const double crystal_size_z, const uint32_t Nyx, 
	const uint32_t Nz, std::vector<double>& store_elements, std::vector<uint32_t>& store_indices, const uint32_t tid, uint32_t& ind, 
	double* elements, mwIndex* indices, const size_t lo, const bool PRECOMPUTE, const double global_factor, uint64_t N2) {

	//const uint32_t zz = z_loop * d_N2 * d_N1;
	const uint32_t apu = perpendicular_start(d_b, d, d_d1, d_N1);
	double jelppi = 0.;
	uint32_t hpk = N2;
	for (int32_t zz = static_cast<int32_t>(z_loop); zz >= 0; zz--) {
		for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
			double local_ele = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
				break;
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
			double local_ele = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
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
	for (uint32_t zz = z_loop + 1u; zz < Nz; zz++) {
		for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
			double local_ele = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
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
			double local_ele = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
			if (local_ele <= THR)
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
	temp *= global_factor;
	if (PRECOMPUTE) {
		for (int32_t zz = static_cast<int32_t>(z_loop); zz >= 0; zz--) {
			for (int32_t uu = static_cast<int32_t>(apu); uu >= 0; uu--) {
				double local_ele = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
				if (local_ele <= THR)
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
				double local_ele = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
				if (local_ele <= THR)
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
				double local_ele = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
				if (local_ele <= THR)
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
				double local_ele = compute_element_orth_3D(detectors, xl, yl, zl, crystal_size_z, center1[uu], center2, z_center[zz]);
				if (local_ele <= THR)
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

// Calculate the normalization factor and RHS in the perpendicular case in orthogonal ray tracer (3D case)
//void orth_distance_rhs_perpendicular_mfree_3D(const double* center1, const double center2, const double* z_center, const double temp, double& ax, 
//	const double d_b, const double d, const double d_d1, const uint32_t d_N1, const uint32_t d_N2, const uint32_t z_loop, const uint32_t d_N, 
//	const uint32_t d_NN, const bool no_norm, double* rhs, double* Summ, const bool RHS, const bool SUMMA, Det detectors, const double xl, const double yl, const double zl, 
//	const double crystal_size_z, const uint32_t Nyx, const uint32_t Nz, const std::vector<double> store_elements, const std::vector<uint32_t> store_indices,
//	const uint32_t tid, uint32_t ind, std::vector<double>& elements, std::vector<uint32_t>& indices, uint32_t N2) {
//
//	//const uint32_t zz = z_loop * d_N2 * d_N1;
//	const uint32_t apu = perpendicular_start(d_b, d, d_d1, d_N1);
//	uint32_t hpk = N2;
//	for (int32_t uu = 0; uu < ind; uu++) {
//		double local_ele = store_elements[tid + uu];
//		uint32_t local_ind = store_indices[tid + uu];
//		local_ele *= temp;
//		if (RHS) {
//			for (uint32_t kk = 0u; kk < d_N2; kk++) {
//#pragma omp atomic
//				rhs[local_ind] += (local_ele * ax);
//				if (no_norm == 0) {
//#pragma omp atomic
//					Summ[local_ind] += local_ele;
//				}
//				local_ind += d_NN;
//			}
//		}
//		else if (SUMMA) {
//			for (uint32_t kk = 0u; kk < d_N2; kk++) {
//#pragma omp atomic
//				Summ[local_ind] += local_ele;
//				local_ind += d_NN;
//			}
//		}
//		else {
//			for (uint32_t kk = 0u; kk < d_N2; kk++) {
//				indices[hpk] = local_ind;
//				elements[hpk] = local_ele;
//				local_ind += d_NN;
//				hpk++;
//			}
//		}
//	}
//}

void setThreads() {
	if (omp_get_max_threads() == 1) {
		int n_threads = std::thread::hardware_concurrency();
		omp_set_num_threads(n_threads);
	}
}