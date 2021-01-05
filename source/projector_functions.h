/**************************************************************************
* Header file for the various functions required by the original Siddon, 
* improved Siddon, orthogonal distance-based and volume-based ray tracers.
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
#pragma once

#ifdef MATLAB
#include "mex.h"
#elif OCTAVE
#include <octave/oct.h>
#endif
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <time.h>
#include <cstdint>
#include <thread>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif

// Normalized distances below this are discarded in orthogonal ray tracer
#define THR 0.01
#define H_THR 0.99
#define TOF_THR 0.0001
#define _2PI 0.3989422804014327
#define TRAPZ_BINS 6.

struct Det {
	double xd, xs, yd, ys, zd, zs;
};

double norm(const double x, const double y, const double z);

double compute_element_orth_3D(Det detectors, const double xl, const double yl, const double zl, const double crystal_size_z,
	const double xp, const double yp, const double zp);

void computeIndices(const bool RHS, const bool SUMMA, const bool OMP, const bool PRECOMP, const bool DISCARD, double local_ele, double& temp, double& ax,
	const bool no_norm, double* Summ, double* rhs, const double local_sino, const double* osem_apu, const uint64_t N2, size_t* indices,
	std::vector<double>& elements, std::vector<uint32_t>& v_indices, size_t& idx, const uint32_t local_ind, const uint64_t N22);

void get_detector_coordinates_raw(const uint32_t det_per_ring, const double* x, const double* y, const double* z, Det& detectors,
	const uint16_t* L, const size_t ll, const uint32_t* pseudos, const uint32_t pRows, const bool list_mode_format = false);

void get_detector_coordinates_raw_N(const uint32_t det_per_ring, const double* x, const double* y, const double* z, Det& detectors,
	const uint16_t* L, const size_t ll, const uint32_t* pseudos, const uint32_t pRows, const uint16_t lor, const double cr_pz, 
	const uint16_t n_rays, const uint16_t n_rays3D);

void get_detector_coordinates_noalloc(const double* x, const double* y, const double* z, const uint32_t size_x, Det& detectors, int& ll, const uint32_t* index,
	int& lz, const uint32_t TotSinos, size_t oo);

void get_detector_coordinates(const double* x, const double* y, const double* z, const uint32_t size_x, Det& detectors, const uint32_t* xy_index,
	const uint16_t* z_index, const uint32_t TotSinos, const size_t oo);

void get_detector_coordinates_mr(const double* x, const double* y, const double* z, const uint32_t size_x, Det& detectors, const uint32_t* xy_index,
	const uint16_t* z_index, const uint32_t TotSinos, const size_t oo, const uint16_t lor, const double cr_pz, const uint16_t n_rays, const uint16_t n_rays3D);

uint32_t z_ring(const double zmax, const double zs, const double NSlices);

void s_g_d(const double tmin, const double t_min, const double tmax, const double t_max, uint32_t& v_min, uint32_t& v_max, double& t_0, int32_t& v_u, 
	const double diff, const double b, const double d, const double s, const uint32_t N);

void d_g_s(const double tmin, const double t_min, const double tmax, const double t_max, uint32_t& v_min, uint32_t& v_max, double& t_0, int32_t& v_u, 
	const double diff, const double b, const double d, const double s, const uint32_t N);

void s_g_d_precomp(const double tmin, const double t_min, const double tmax, const double t_max, uint32_t& v_min, uint32_t& v_max, double& t_0, int32_t& v_u, 
	const double diff, const double b, const double d, const double s, const uint32_t N);

void d_g_s_precomp(const double tmin, const double t_min, const double tmax, const double t_max, uint32_t& v_min, uint32_t& v_max, double& t_0, int32_t& v_u, 
	const double diff, const double b, const double d, const double s, const uint32_t N);

double pixel_value(const double t, const double tc, const double L);

void compute_attenuation(double& tc, double& jelppi, const double LL, const double t0, const int tempi, const int tempj, const int tempk, const uint32_t Nx, 
	const uint32_t Nyx, const double* atten);

void att_corr_vec(const std::vector<double> templ_ijk, const std::vector<uint32_t> temp_koko, const double* atten, double& temp, const size_t Np);

void att_corr_vec_precomp(const double* elements, const double* atten, const size_t* indices, const size_t Np, const uint64_t N2, double& temp);

double att_corr_scalar(double templ_ijk, uint32_t tempk, const double* atten, double& temp, const uint32_t N1, const uint32_t N);

void att_corr_scalar_orth(uint32_t tempk, const double* atten, double& temp, const uint32_t N1, const uint32_t N2, const double d);

double perpendicular_elements(const uint32_t N, const double dd, const std::vector<double> vec, const double d, const uint32_t z_ring, const uint32_t N1, 
	const uint32_t N2, const double* atten, const double* norm_coef, const bool attenuation_correction, const bool normalization, uint32_t& tempk, 
	const uint32_t NN, const size_t lo, const double global_factor, const bool scatter, const double* scatter_coef);

double perpendicular_elements_multiray(const uint32_t N, const double dd, const std::vector<double> vec, const double d, const uint32_t z_ring,
	const uint32_t N1, const uint32_t N2, const double* atten, const bool attenuation_correction, int32_t& tempk, const uint32_t NN, double& jelppi);

// this function was taken from: https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
// Currently unused
template<typename T>
inline std::vector<size_t> sort_indexes(const std::vector<T>& v)
{
	// initialize original index locations
	std::vector<size_t> idx(v.size());
	std::iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

	return idx;
}

int32_t voxel_index(const double pt, const double diff, const double d, const double apu);

void orth_perpendicular(const uint32_t N, const double dd, const std::vector<double> vec, const double d, const uint32_t z_ring, const uint32_t N1, 
	const uint32_t N2, const double diff2, const double* center1, const double crystal_size, const double length, std::vector<double>& vec1, 
	std::vector<double>& vec2, double& temp, uint32_t& tempk); 

void orth_perpendicular_3D(const double dd, const std::vector<double> vec, const uint32_t z_ring, const uint32_t N1, const uint32_t N2, const uint32_t Nz, 
	const uint32_t Nyx, const uint32_t d_N, const uint32_t d_NN, const Det detectors, const double xl, const double yl, const double zl, const double* center1, 
	const double center2, const double* z_center, const double crystal_size_z, int& hpk, double& temp, uint32_t& tempk, size_t* indices, double* elements, 
	const uint64_t Np);

void orth_perpendicular_np_3D(const double dd, const std::vector<double> vec, const uint32_t z_ring, const uint32_t N1, const uint32_t N2, const uint32_t Nz, 
	const uint32_t Nyx, const uint32_t d_N, const uint32_t d_NN, const Det detectors, const double xl, const double yl, const double zl, const double* center1, 
	const double center2, const double* z_center, const double crystal_size_z, int& hpk, double& temp, uint32_t& tempk, std::vector<uint32_t>& indices, 
	std::vector<double>& elements);

void orth_perpendicular_precompute(const uint32_t N1, const uint32_t N2, const double dd, const std::vector<double> vec, 
	const double* center1, const double center2, const double* z_center, const double kerroin, size_t& temp_koko, const Det detectors,
	const double xl, const double yl, const double zl, const int32_t tempk);

void orth_perpendicular_precompute_3D(const uint32_t N1, const uint32_t N2, const uint32_t Nz, const double dd, const std::vector<double> vec,
	const double* center1, const double center2, const double* z_center, const double crystal_size_z, size_t& temp_koko, const Det detectors,
	const double xl, const double yl, const double zl, const uint32_t z_loop);

void orth_distance_precompute(const int32_t tempi, const uint32_t Nx, const double y_diff, const double x_diff, const double* x_center, const double y_center, 
	const double kerroin, const double length, uint16_t& temp_koko);

void orth_distance_precompute_3D(const int32_t tempi, const uint32_t N, const uint32_t Nz, const double y_diff, const double x_diff, const double z_diff,
	const double* x_center, const double y_center, const double* z_center, const double crystal_size_z, const Det detectors,
	uint16_t& temp_koko, const int32_t tempk, const size_t lo, const int32_t n_tempk, const int32_t dec, const int32_t decx);

void orth_distance_full(const int32_t tempi, const uint32_t Nx, const double y_diff, const double x_diff, const double z_diff, const double* y_center, const double* x_center, const double* z_center,
	const Det detectors, const double kerroin, double& temp, const uint32_t tempijk, const uint32_t NN, const int32_t tempj, const int32_t tempk,
	const double local_sino, double& ax, const double* osem_apu, const bool no_norm, const bool RHS, const bool SUMMA, const bool OMP, const bool PRECOMP,
	const bool DISCARD, double* rhs, double* Summ, size_t* indices, std::vector<double>& elements, std::vector<uint32_t>& v_indices, size_t& idx, const uint32_t Ny,
	const uint32_t N1, const int8_t start, int32_t ju, uint8_t xyz, std::vector<double>& store_elements, std::vector<uint32_t>& store_indices,
	const uint32_t dec_v, uint32_t& ind, uint64_t N2 = 0ULL, uint64_t N22 = 0ULL);

void orth_distance_3D_full(int32_t tempi, const uint32_t Nx, const uint32_t Nz, const double y_diff, const double x_diff, const double z_diff,
	const double* y_center, const double* x_center, const double* z_center, double& temp, const uint32_t NN, int32_t tempj, int32_t tempk,
	const double local_sino, double& ax, const double* osem_apu, const Det detectors, const uint32_t Nyx, const double kerroin,
	const bool no_norm, const bool RHS, const bool SUMMA, const bool OMP, const bool PRECOMP, const bool DISCARD, double* rhs, double* Summ, size_t* indices,
	std::vector<double>& elements, std::vector<uint32_t>& v_indices, size_t& idx, const uint32_t Ny, const uint32_t N1, const int start,
	const int32_t iu, const int32_t ju, const int loppu, std::vector<double>& store_elements, std::vector<uint32_t>& store_indices,
	const uint32_t dec_v, uint32_t& ind, uint64_t N2 = 0ULL, uint64_t N22 = 0ULL);

void volume_distance_3D_full(int32_t tempi, const uint32_t Nx, const uint32_t Nz, const double y_diff, const double x_diff, const double z_diff,
	const double* y_center, const double* x_center, const double* z_center, double& temp, const uint32_t NN, int32_t tempj, int32_t tempk,
	const double local_sino, double& ax, const double* osem_apu, const Det detectors, const uint32_t Nyx, const double kerroin,
	const bool no_norm, const bool RHS, const bool SUMMA, const bool OMP, const bool PRECOMP, const bool DISCARD, double* rhs, double* Summ, size_t* indices,
	std::vector<double>& elements, std::vector<uint32_t>& v_indices, size_t& idx, const uint32_t Ny, const uint32_t N1, const int start,
	const int32_t iu, const int32_t ju, const int loppu, std::vector<double>& store_elements, std::vector<uint32_t>& store_indices,
	const uint32_t tid, uint32_t& ind, const double bmax, const double bmin, const double Vmax, const double* V, uint64_t N2 = 0ULL, uint64_t N22 = 0ULL);

void volume_perpendicular_precompute(const uint32_t N1, const uint32_t N2, const uint32_t Nz, const double dd, const std::vector<double> vec,
	const double* center1, const double center2, const double* z_center, const double crystal_size_z, size_t& temp_koko, const Det detectors,
	const double xl, const double yl, const double zl, const uint32_t z_loop, const double bmax, const double bmin, const double Vmax,
	const double* V);

void nominator_mfree(double& ax, const double Sino, const double epps, const double temp, const bool randoms_correction, const double* randoms, 
	const size_t lo);

double compute_element_orth_mfree(const double x_diff, const double y_diff, const double x_center, const double length_);

uint32_t compute_ind_orth_mfree(const uint32_t tempi, const uint32_t temp_ijk, const uint32_t d_N);

uint32_t compute_ind_orth_mfree_3D(const uint32_t tempi, const uint32_t tempijk, const uint32_t tempk, const uint32_t d_N, const uint32_t Nyx);

void denominator_mfree(const double local_ele, double& axOSEM, const double d_OSEM);

uint32_t perpendicular_start(const double d_b, const double d, const double d_d, const uint32_t d_N);

void orth_distance_denominator_perpendicular_mfree(const double* center1, const double center2, const double* z_center, const double kerroin,
	double& temp, const bool d_attenuation_correction, const bool normalization, double& ax, const double d_b, const double d, const double d_d1,
	const uint32_t d_N1, const uint32_t d_N2, const uint32_t z_loop, const double* d_atten, const double* norm_coef, const double local_sino, const uint32_t d_N, const uint32_t d_NN,
	const double* d_OSEM, const Det detectors, const double xl, const double yl, const double zl, std::vector<double>& store_elements, std::vector<uint32_t>& store_indices,
	const uint32_t dec_v, uint32_t& ind, double* elements, size_t* indices, const size_t lo, const bool PRECOMPUTE, const double global_factor, 
	const bool scatter, const double* scatter_coef, const uint64_t N2 = 0ULL);

void volume_distance_denominator_perpendicular_mfree_3D(const double* center1, const double center2, const double* z_center, double& temp,
	const bool d_attenuation_correction, const bool normalization, double& ax, const double d_b, const double d, const double d_d1, const uint32_t d_N1,
	const uint32_t d_N2, const uint32_t z_loop, const double* d_atten, const double* norm_coef, const double local_sino, const uint32_t d_N, const uint32_t d_NN,
	const double* d_OSEM, Det detectors, const double xl, const double yl, const double zl, const double crystal_size_z, const uint32_t Nyx,
	const uint32_t Nz, std::vector<double>& store_elements, std::vector<uint32_t>& store_indices, const uint32_t tid, uint32_t& ind,
	double* elements, size_t* indices, const size_t lo, const bool PRECOMPUTE, const double global_factor, const double bmax, const double bmin, const double Vmax,
	const double* V, const bool scatter, const double* scatter_coef, const uint64_t N2 = 0ULL);

void orth_distance_rhs_perpendicular_mfree(const double* center1, const double center2, const double* z_center, const double kerroin,
	const double temp, double& ax, const double d_b, const double d, const double d_d1, const uint32_t d_N1, const uint32_t d_N2,
	const uint32_t z_loop, const uint32_t d_N, const uint32_t d_NN, const bool no_norm, double* rhs, double* Summ, const bool RHS, const bool SUMMA,
	const Det detectors, const double xl, const double yl, const double zl, const std::vector<double> store_elements, const std::vector<uint32_t> store_indices,
	const uint32_t tid, uint32_t ind, double* elements, size_t* indices, uint64_t N2 = 0ULL);

//void orth_distance_summ_perpendicular_mfree(const double diff2, const double* center1, const double kerroin, const double length_, const double temp,
//	double ax, const double d_b, const double d, const double d_d1, const uint32_t d_N1, const uint32_t d_N2, const uint32_t z_loop, const uint32_t d_N, 
//	const uint32_t d_NN, double* Summ);

//int orth_siddon_no_precompute(const uint32_t loop_var_par, const uint32_t size_x, const double zmax, const uint32_t TotSinos, std::vector<uint32_t>& indices, 
//	std::vector<double>& elements, uint16_t* lor, const double maxyy, const double maxxx, const std::vector<double>& xx_vec, const double dy,
//	const std::vector<double>& yy_vec, const double* atten, const double* norm_coef, const double* x, const double* y, const double* z_det, 
//	const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const double dx, const double dz, const double bx, const double by, 
//	const double bz, const uint32_t* index, const bool attenuation_correction, const bool normalization, const bool raw, const uint32_t det_per_ring, 
//	const uint32_t blocks, const uint32_t block1, const uint16_t* L, const uint32_t* pseudos, const uint32_t pRows, const double crystal_size, 
//	const double crystal_size_z, const double* y_center, const double* x_center, const double* z_center, const uint32_t dec_v);

void orth_siddon_precomputed(const int64_t loop_var_par, const uint32_t size_x, const double zmax, size_t* indices, double* elements, const double maxyy,
	const double maxxx, const std::vector<double>& xx_vec, const double dy, const std::vector<double>& yy_vec, const double* atten, const double* norm_coef,
	const double* x, const double* y, const double* z_det, const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const double dx, 
	const double dz, const double bx, const double by, const double bz, const bool attenuation_correction, const bool normalization, const uint16_t* lor1, 
	const uint64_t* lor2, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t TotSinos, const uint16_t* L, const uint32_t* pseudos, 
	const uint32_t pRows, const uint32_t det_per_ring, const bool raw, const bool attenuation_phase, double* length, const double crystal_size, 
	const double crystal_size_z, double* y_center, double* x_center, const double* z_center, const double global_factor, const bool scatter, 
	const double* scatter_coef, const uint32_t nCores = 1U);

void vol_siddon_precomputed(const int64_t loop_var_par, const uint32_t size_x, const double zmax, size_t* indices, double* rhs, const double maxyy,
	const double maxxx, const std::vector<double>& xx_vec, const double dy, const std::vector<double>& yy_vec, const double* atten, const double* norm_coef,
	const double* x, const double* y, const double* z_det, const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const double dx,
	const double dz, const double bx, const double by, const double bz, const bool attenuation_correction, const bool normalization, const uint16_t* lor1,
	const uint64_t* lor2, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t TotSinos, const uint16_t* L, const uint32_t* pseudos,
	const uint32_t pRows, const uint32_t det_per_ring, const bool raw, const bool attenuation_phase, double* length, const double crystal_size,
	const double crystal_size_z, double* y_center, double* x_center, const double* z_center, const double global_factor, const double bmin,
	const double bmax, const double Vmax, const double* V, const bool scatter, const double* scatter_coef, const uint32_t nCores = 1U);

int improved_siddon_no_precompute(const int64_t loop_var_par, const uint32_t size_x, const double zmax, const uint32_t TotSinos, std::vector<uint32_t>& indices,
	std::vector<double>& elements, uint16_t* lor, const double maxyy, const double maxxx, const std::vector<double>& xx_vec, const double dy,
	const std::vector<double>& yy_vec, const double* atten, const double* norm_coef, const double* x, const double* y, const double* z_det, 
	const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const double dx, const double dz, const double bx, const double by, 
	const double bz, const uint32_t* index, const bool attenuation_correction, const bool normalization, const bool raw, const uint32_t det_per_ring, 
	const uint32_t blocks, const uint32_t block1, const uint16_t* L, const uint32_t* pseudos, const uint32_t pRows, const double global_factor, 
	const bool scatter, const double* scatter_coef);

int original_siddon_no_precompute(const int64_t loop_var_par, const uint32_t size_x, const double zmax, const uint32_t TotSinos, std::vector<uint32_t>& indices,
	std::vector<double>& elements, uint16_t* lor, const double maxyy, const double maxxx, const std::vector<double>& xx_vec, const double dy,
	const std::vector<double>& yy_vec, const double* atten, const double* norm_coef, const double* x, const double* y, const double* z_det, 
	const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const double dx, const double dz, const double bx, const double by, 
	const double bz, const uint32_t* index, const bool attenuation_correction, const bool normalization, const bool raw, const uint32_t det_per_ring, 
	const uint32_t blocks, const uint32_t block1, const uint16_t* L, const uint32_t* pseudos, const uint32_t pRows, const std::vector<double>& iij_vec, 
	const std::vector<double>& jjk_vec, const std::vector<double>& kkj_vec, const double global_factor, const bool scatter, const double* scatter_coef);

extern void improved_siddon_precomputed(const int64_t loop_var_par, const uint32_t size_x, const double zmax, size_t* indices, double* elements, const double maxyy,
	const double maxxx, const std::vector<double>& xx_vec, const double dy, const std::vector<double>& yy_vec, const double* atten, const double* norm_coef, 
	const double* x, const double* y, const double* z_det, const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const double dx, 
	const double dz, const double bx, const double by, const double bz, const bool attenuation_correction, const bool normalization, const uint16_t* lor1, 
	const uint64_t* lor2, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t TotSinos, const uint16_t* L, const uint32_t* pseudos, 
	const uint32_t pRows, const uint32_t det_per_ring, const bool raw, const bool attenuation_phase, double* length, const double global_factor, 
	const bool scatter, const double* scatter_coef, const uint32_t nCores = 1U);

void sequential_improved_siddon(const int64_t loop_var_par, const uint32_t size_x, const double zmax, double* Summ, double* rhs, const double maxyy,
	const double maxxx, const std::vector<double>& xx_vec, const double dy, const std::vector<double>& yy_vec, const double* atten, const double* norm_coef, 
	const double* randoms, const double* x, const double* y, const double* z_det, const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny, 
	const uint32_t Nz, const double dx, const double dz, const double bx, const double by, const double bz,	const bool attenuation_correction, 
	const bool normalization, const bool randoms_correction, const uint16_t* lor1, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t TotSinos,
	const double epps, const double* Sino, double* osem_apu, const uint16_t* L, const uint32_t* pseudos, const uint32_t pRows, const uint32_t det_per_ring,
	const bool raw, const bool no_norm, const double global_factor, const uint8_t fp, const bool scatter, const double* scatter_coef, const bool TOF,
	const int64_t TOFSize, const double sigma_x, const double* TOFCenter, const int64_t nBins, const uint32_t dec_v, const uint32_t nCores = 1U);

void sequential_improved_siddon_no_precompute(const int64_t loop_var_par, const uint32_t size_x, const double zmax, double* Summ, double* rhs, const double maxyy,
	const double maxxx, const std::vector<double>& xx_vec, const double dy, const std::vector<double>& yy_vec, const double* atten, const double* norm_coef, 
	const double* randoms, const double* x, const double* y, const double* z_det, const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny, 
	const uint32_t Nz, const double dx, const double dz, const double bx, const double by, const double bz, const bool attenuation_correction, 
	const bool normalization, const bool randoms_correction, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t TotSinos, const double epps, 
	const double* Sino, double* osem_apu, const uint16_t* L, const uint32_t* pseudos, const size_t pRows, const uint32_t det_per_ring, const bool raw, 
	const double cr_pz, const bool no_norm, const uint16_t n_rays, const uint16_t n_rays3D, const double global_factor, const uint8_t fp, const bool list_mode_format,
	const bool scatter, const double* scatter_coef, const bool TOF, const int64_t TOFSize, const double sigma_x, const double* TOFCenter,
	const int64_t nBins, const uint32_t dec_v, const uint32_t nCores = 1U);

void sequential_orth_siddon(const int64_t loop_var_par, const uint32_t size_x, const double zmax, double* Summ, double* rhs, const double maxyy,
	const double maxxx,	const std::vector<double>& xx_vec, const double dy, const std::vector<double>& yy_vec, const double* atten, const double* norm_coef, 
	const double* randoms, const double* x, const double* y, const double* z_det, const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, 
	const double dx, const double dz, const double bx, const double by, const double bz, const bool attenuation_correction, const bool normalization, 
	const bool randoms_correction, const uint16_t* lor1, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t TotSinos, const double epps, 
	const double* Sino, double* osem_apu, const uint16_t* L, const uint32_t* pseudos, const uint32_t pRows, const uint32_t det_per_ring, const bool raw, 
	const double crystal_size_xy, double* x_center, double* y_center, const double* z_center, const double crystal_size_z, const bool no_norm, 
	const uint32_t dec_v, const double global_factor, const uint8_t fp, const bool scatter, const double* scatter_coef, const uint32_t nCores = 1U);

void sequential_orth_siddon_no_precomp(const int64_t loop_var_par, const uint32_t size_x, const double zmax, double* Summ, double* rhs, const double maxyy,
	const double maxxx, const std::vector<double>& xx_vec, const double dy, const std::vector<double>& yy_vec, const double* atten, const double* norm_coef, 
	const double* randoms, const double* x, const double* y, const double* z_det, const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny, 
	const uint32_t Nz, const double dx, const double dz, const double bx, const double by, const double bz, const bool attenuation_correction, 
	const bool normalization, const bool randoms_correction, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t TotSinos, const double epps, 
	const double* Sino, double* osem_apu, const uint16_t* L, const uint32_t* pseudos, const uint32_t pRows, const uint32_t det_per_ring, const bool raw, 
	const double crystal_size_xy, double* x_center, double* y_center, const double* z_center, const double crystal_size_z, const bool no_norm, 
	const uint32_t dec_v, const double global_factor, const uint8_t fp, const bool list_mode_format, const bool scatter, const double* scatter_coef, 
	const uint32_t nCores = 1U);


// Source: http://www.alecjacobson.com/weblog/?p=4544 & https://ideone.com/Z7zldb
// Unused, replaced with OpenMP
//class ThreadPool {
//
//public:
//
//	template<typename Index, typename Callable>
//	static void ParallelFor(Index start, Index end, Callable func) {
//		// Estimate number of threads in the pool
//		const static unsigned nb_threads_hint = std::thread::hardware_concurrency();
//		const static unsigned nb_threads = (nb_threads_hint == 0u ? 8u : nb_threads_hint);
//
//		// Size of a slice for the range functions
//		Index n = end - start + 1;
//		Index slice = (Index)std::round(n / static_cast<double> (nb_threads));
//		slice = std::max(slice, Index(1));
//
//		// [Helper] Inner loop
//		auto launchRange = [&func](uint32_t k1, uint32_t k2) {
//			for (Index k = k1; k < k2; k++) {
//				func(k);
//			}
//		};
//
//		// Create pool and launch jobs
//		std::vector<std::thread> pool;
//		pool.reserve(nb_threads);
//		Index i1 = start;
//		Index i2 = std::min(start + slice, end);
//		for (unsigned i = 0; i + 1 < nb_threads && i1 < end; ++i) {
//			pool.emplace_back(launchRange, i1, i2);
//			i1 = i2;
//			i2 = std::min(i2 + slice, end);
//		}
//		if (i1 < end) {
//			pool.emplace_back(launchRange, i1, end);
//		}
//
//		// Wait for jobs to finish
//		for (std::thread &t : pool) {
//			if (t.joinable()) {
//				t.join();
//			}
//		}
//	}
//
//	// Serial version for easy comparison
//	template<typename Index, typename Callable>
//	static void SequentialFor(Index start, Index end, Callable func) {
//		for (Index i = start; i < end; i++) {
//			func(i);
//		}
//	}
//
//};

bool siddon_pre_loop_2D(const double b1, const double b2, const double diff1, const double diff2, const double max1, const double max2,
	const double d1, const double d2, const uint32_t N1, const uint32_t N2, int32_t& temp1, int32_t& temp2, double& t1u, double& t2u, uint32_t& Np,
	const int TYPE, const double ys, const double xs, const double yd, const double xd, double& tc, int32_t& u1, int32_t& u2, double& t10, double& t20);

bool siddon_pre_loop_3D(const double bx, const double by, const double bz, const double x_diff, const double y_diff, const double z_diff,
	const double maxxx, const double maxyy, const double bzb, const double dx, const double dy, const double dz,
	const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, int32_t& tempi, int32_t& tempj, int32_t& tempk, double& tyu, double& txu, double& tzu,
	uint32_t& Np, const int TYPE, const Det detectors, double& tc, int32_t& iu, int32_t& ju, int32_t& ku, double& tx0, double& ty0, double& tz0);

void orth_distance_denominator_perpendicular_mfree_3D(const double* center1, const double center2, const double* z_center, double& temp,
	const bool d_attenuation_correction, const bool normalization, double& ax, const double d_b, const double d, const double d_d1, const uint32_t d_N1,
	const uint32_t d_N2, const uint32_t z_loop, const double* d_atten, const double* norm_coef, const double local_sino, const uint32_t d_N, const uint32_t d_NN,
	const double* d_OSEM, Det detectors, const double xl, const double yl, const double zl, const double crystal_size_z, const uint32_t Nyx,
	const uint32_t Nz, std::vector<double>& store_elements, std::vector<uint32_t>& store_indices, const uint32_t tid, uint32_t& ind, 
	double* elements, size_t* indices, const size_t lo, const bool PRECOMPUTE, const double global_factor, const bool scatter, const double* scatter_coef, 
	uint64_t N2 = 0ULL);

void improved_siddon_precomputation_phase(const int64_t loop_var_par, const uint32_t size_x, const double zmax, const uint32_t TotSinos, uint16_t* lor, 
	const double maxyy,	const double maxxx, const std::vector<double>& xx_vec, const std::vector<double>& z_det_vec, const double dy, const std::vector<double>& yy_vec,
	const double* x, const double* y, const double* z_det, const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const double dx, 
	const double dz,	const double bx, const double by, const double bz, const uint32_t block1, const uint32_t blocks, const uint16_t* L, const uint32_t* pseudos,
	const bool raw, const uint32_t pRows, const uint32_t det_per_ring, const uint32_t type, uint16_t* lor_orth, uint16_t* lor_vol, const double crystal_size, 
	const double crystal_size_z,	const double* x_center, const double* y_center, const double* z_center, const double bmin, const double bmax, const double Vmax, 
	const double* V, const uint32_t nCores = 1U);

void sequential_volume_siddon_no_precomp(const int64_t loop_var_par, const uint32_t size_x, const double zmax, double* Summ, double* rhs, const double maxyy,
	const double maxxx, const std::vector<double>& xx_vec, const double dy, const std::vector<double>& yy_vec, const double* atten, const double* norm_coef,
	const double* randoms, const double* x, const double* y, const double* z_det, const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny,
	const uint32_t Nz, const double dx, const double dz, const double bx, const double by, const double bz, const bool attenuation_correction,
	const bool normalization, const bool randoms_correction, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t TotSinos,
	const double epps, const double* Sino, double* osem_apu, const uint16_t* L, const uint32_t* pseudos, const uint32_t pRows, const uint32_t det_per_ring,
	const bool raw, const double Vmax, double* x_center, double* y_center, const double* z_center, const double bmin, const double bmax, const double* V,
	const bool no_norm, const uint32_t dec_v, const double global_factor, const uint8_t fp, const bool list_mode_format, const bool scatter, const double* scatter_coef, 
	const uint32_t nCores = 1U);

void sequential_volume_siddon(const int64_t loop_var_par, const uint32_t size_x, const double zmax, double* Summ, double* rhs, const double maxyy,
	const double maxxx, const std::vector<double>& xx_vec, const double dy, const std::vector<double>& yy_vec, const double* atten, const double* norm_coef,
	const double* randoms, const double* x, const double* y, const double* z_det, const uint32_t NSlices, const uint32_t Nx, const uint32_t Ny,
	const uint32_t Nz, const double dx, const double dz, const double bx, const double by, const double bz, const bool attenuation_correction,
	const bool normalization, const bool randoms_correction, const uint16_t* lor1, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t TotSinos,
	const double epps, const double* Sino, double* osem_apu, const uint16_t* L, const uint32_t* pseudos, const uint32_t pRows, const uint32_t det_per_ring,
	const bool raw, const double Vmax, double* x_center, double* y_center, const double* z_center, const double bmin, const double bmax, const double* V,
	const bool no_norm, const uint32_t dec_v, const double global_factor, const uint8_t fp, const bool scatter, const double* scatter_coef, const uint32_t nCores = 1U);

//void orth_distance_rhs_perpendicular_mfree_3D(const double* center1, const double center2, const double* z_center, const double temp, double& ax,
//	const double d_b, const double d, const double d_d1, const uint32_t d_N1, const uint32_t d_N2, const uint32_t z_loop, const uint32_t d_N,
//	const uint32_t d_NN, const bool no_norm, double* rhs, double* Summ, const bool RHS, Det detectors, const double xl, const double yl, const double zl,
//	const double crystal_size_z, const uint32_t Nyx, const uint32_t Nz, const std::vector<double> store_elements, const std::vector<uint32_t> store_indices,
//	const uint32_t tid, uint32_t ind);
//
//void orth_distance_summ_perpendicular_mfree_3D(const double* center1, const double center2, const double* z_center, const double temp,
//	double ax, const double d_b, const double d, const double d_d1, const uint32_t d_N1, const uint32_t d_N2, const uint32_t z_loop, const uint32_t d_N, 
//	const uint32_t d_NN, double* Summ, Det detectors, const double xl, const double yl, const double zl, const double crystal_size_z, const uint32_t Nyx, 
//	const uint32_t Nz);

void setThreads();

template <typename T>
T normPDF(const T x, const T mu, const T sigma) {

	const T a = (x - mu) / sigma;

	return (_2PI / sigma * std::exp(-0.5 * a * a));
}

template <typename T>
void TOFLoop(T& TOFSum, const T DD, const int64_t nBins, const T element, std::vector<T>& TOFVal, const T* TOFCenter, 
	const T sigma_x, T& D, const int64_t tid, const T epps) {
	const T dX = element / static_cast<T>(TRAPZ_BINS);
	if (DD > 0) {
		for (int64_t to = 0LL; to < nBins; to++) {
			//TOFVal[to + tid] = (dX * ((normPDF(D, TOFCenter[to], sigma_x) + normPDF(D - element, TOFCenter[to], sigma_x)) / 2.));
			TOFVal[to + tid] = normPDF(D, TOFCenter[to], sigma_x);
			for (int64_t tr = 1LL; tr < static_cast<int64_t>(TRAPZ_BINS) - 1; tr++)
				TOFVal[to + tid] += (normPDF(D - dX * static_cast<T>(tr), TOFCenter[to], sigma_x) * 2.);
			TOFVal[to + tid] += normPDF(D - element, TOFCenter[to], sigma_x);
			TOFVal[to + tid] *= dX;
			TOFSum += TOFVal[to + tid];
		}
		D -= element;
	}
	else {
		for (int64_t to = 0LL; to < nBins; to++) {
			//TOFVal[to + tid] = (dX * ((normPDF(D, TOFCenter[to], sigma_x) + normPDF(D + element, TOFCenter[to], sigma_x)) / 2.));
			TOFVal[to + tid] = normPDF(D, TOFCenter[to], sigma_x);
			for (int64_t tr = 1LL; tr < TRAPZ_BINS - 1; tr++)
				TOFVal[to + tid] += (normPDF(D + dX * static_cast<T>(tr), TOFCenter[to], sigma_x) * 2.);
			TOFVal[to + tid] += normPDF(D + element, TOFCenter[to], sigma_x);
			TOFVal[to + tid] *= dX;
			TOFSum += TOFVal[to + tid];
		}
		D += element;
	}
}

template <typename T>
void TOFWeightsFP(const T DD, const int64_t nBins, const T element, std::vector<T>& TOFVal, const T* TOFCenter, const T sigma_x,
	T& D, const T* osem_apu, const uint32_t tempijk, std::vector<T>& ax, const T epps, const int64_t tid) {
	T TOFSum = 0.;
	TOFLoop(TOFSum, DD, nBins, element, TOFVal, TOFCenter, sigma_x, D, tid, epps);
	T apu = element * osem_apu[tempijk];
	if (TOFSum < epps)
		TOFSum = epps;
	for (int64_t to = 0LL; to < nBins; to++) {
		TOFVal[to + tid] /= TOFSum;
		ax[to] += (apu * TOFVal[to + tid]);
		//ax[to] += (apu * (TOFVal[to] / TOFSum));
	}
}

template <typename T>
void ForwardProject(T& t0, T& tc, const T tu, const T LL, const bool attenuation_correction, T& jelppi, const T* atten, 
	uint32_t& tempijk, const bool TOF, const T DD, const int64_t nBins, std::vector<T>& TOFVal, const T* TOFCenter, 
	const T sigma_x, T& D, const T* osem_apu, std::vector<T>& ax, const T epps, T& temp, int32_t& tempInd, const int32_t u, 
	const uint32_t incr, const int64_t tid, const uint32_t ind, const uint8_t fp = 0) {

	T element = (t0 - tc) * LL;

	if (attenuation_correction)
		jelppi += (element * -atten[tempijk]);

	if (fp != 2) {
		if (TOF) {
			TOFWeightsFP(DD, nBins, element, TOFVal, TOFCenter, sigma_x, D, osem_apu, tempijk, ax, epps, tid + ind * nBins);
		}
		else
			ax[0] += (element * osem_apu[tempijk]);
	}

	temp += element;

	tempInd += u;
	if (u > 0)
		tempijk += incr;
	else
		tempijk -= incr;

	tc = t0;
	t0 += tu;
}

template <typename T>
void TOFWeightsSumm(const T DD, const int64_t nBins, T element, std::vector<T>& TOFVal, const T* TOFCenter, const T sigma_x,
	T& D, const T epps, const T temp, T& val, const int64_t tid) {
	//T TOFSum = epps;
	//TOFLoop(TOFSum, DD, nBins, element, TOFVal, TOFCenter, sigma_x, D, tid);
	element *= temp;
	val = 0.;
	for (int64_t to = 0LL; to < nBins; to++) {
		//val += element * (TOFVal[to] / TOFSum);
		val += element * TOFVal[to + tid];
	}
}

template <typename T>
T TOFWeightsBP(const T DD, const int64_t nBins, T element, std::vector<T>& TOFVal, const T* TOFCenter, const T sigma_x,
	T& D, const std::vector<T>& yax, const T epps, const T temp, T& val, T* rhs, const int64_t tid) {
	//T TOFSum = epps;
	//TOFLoop(TOFSum, DD, nBins, element, TOFVal, TOFCenter, sigma_x, D, tid);
	element *= temp;
	val = 0.;
	T yaxTOF = 0.;
	for (int64_t to = 0LL; to < nBins; to++) {
		//T apu = element * (TOFVal[to] / TOFSum);
		T apu = element * TOFVal[to + tid];
		val += apu;
		yaxTOF += apu * yax[to];
	}
	return yaxTOF;
}

template <typename T>
void backwardProjection(T& t0, T& tc, const T tu, const T LL, uint32_t& tempijk, const bool TOF, const T DD, const int64_t nBins, 
	std::vector<T>& TOFVal, const T* TOFCenter, const T sigma_x, T& D, std::vector<T>& yax, const T epps, T& temp, 
	const int32_t u, const uint32_t incr, const bool no_norm, T* rhs, T* Summ, const int64_t tid, const uint32_t ind) {
	T element = (t0 - tc) * LL;
	T val = element;
	T val_rhs = 0.;

	if (TOF) {
		val_rhs = TOFWeightsBP(DD, nBins, element, TOFVal, TOFCenter, sigma_x, D, yax, epps, temp, val, rhs, tid + ind * nBins);
	}
	else {
		val *= temp;
		val_rhs = val * yax[0];
	}
#pragma omp atomic
	rhs[tempijk] += (val_rhs);
	if (no_norm == 0 && val > 0.) {
#pragma omp atomic
		Summ[tempijk] += val;
	}

	if (u > 0)
		tempijk += incr;
	else
		tempijk -= incr;

	tc = t0;
	t0 += tu;
}

template <typename T>
void sensImage(T& t0, T& tc, const T tu, const T LL, uint32_t& tempijk, const bool TOF, const T DD, const int64_t nBins,
	std::vector<T>& TOFVal, const T* TOFCenter, const T sigma_x, T& D, const T epps, T& temp, 
	const int32_t u, const uint32_t incr, const bool no_norm, T* Summ, const int64_t tid, const uint32_t ind) {
	T element = (t0 - tc) * LL;
	T val = element;

	if (TOF) {
		TOFWeightsSumm(DD, nBins, element, TOFVal, TOFCenter, sigma_x, D, epps, temp, val, tid + ind * nBins);
	}
	else
		val *= temp;

	if (no_norm == 0 && val > 0.) {
#pragma omp atomic
		Summ[tempijk] += val;
	}

	if (u > 0)
		tempijk += incr;
	else
		tempijk -= incr;

	tc = t0;
	t0 += tu;
}

template <typename T>
void TOFDis(const T x_diff, const T y_diff, const T z_diff, const T tc, const T LL, T& D, T& DD) {
	const T xI = x_diff * tc;
	const T yI = y_diff * tc;
	const T zI = z_diff * tc;
	D = std::sqrt(xI * xI + yI * yI + zI * zI) - LL / 2.;
	DD = D;
}