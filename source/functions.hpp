/**************************************************************************
* Header for ArrayFire functions. Used by implementation 2.
*
* Copyright(C) 2020 Ville - Veikko Wettenhovi
*
* This program is free software : you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <https://www.gnu.org/licenses/>.
***************************************************************************/
#pragma once
#include <arrayfire.h>
#include <algorithm>
#include <vector>
#include <cmath>
#ifdef OPENCL
#include "precomp.h"
#include <af/opencl.h>
#else
#include <nvrtc.h>
#include <cublas.h>
#include <cuda.h>
#include <af/cuda.h>
#include <mex.h>
#endif
#define DEBUG false

#pragma pack(1)

// Struct for the TV-prior
typedef struct TVdata_ {
	af::array s1, s2, s3, s4, s5, s6, s7, s8, s9, reference_image, APLSReference;
	bool TV_use_anatomical;
	float tau, TVsmoothing, T, C, eta, APLSsmoothing, TGVAlpha, TGVBeta, SATVPhi = 0.f;
	uint32_t TVtype = 0;
	uint32_t NiterTGV;
} TVdata;

// Struct for the various estimates
// Default values are scalars (needed for OpenCL kernel)
typedef struct AF_im_vectors_ {
	af::array OSEM, MLEM, RAMLA, MRAMLA, ROSEM, RBI, DRAMA, COSEM, ECOSEM, ACOSEM, 
		MRP_OSEM, MRP_MLEM, MRP_MBSREM, MRP_BSREM, MRP_ROSEM, MRP_RBI, MRP_COSEM,
		Quad_OSEM, Quad_MLEM, Quad_MBSREM, Quad_BSREM, Quad_ROSEM, Quad_RBI, Quad_COSEM,
		Huber_OSEM, Huber_MLEM, Huber_MBSREM, Huber_BSREM, Huber_ROSEM, Huber_RBI, Huber_COSEM,
		L_OSEM, L_MLEM, L_MBSREM, L_BSREM, L_ROSEM, L_RBI, L_COSEM,
		FMH_OSEM, FMH_MLEM, FMH_MBSREM, FMH_BSREM, FMH_ROSEM, FMH_RBI, FMH_COSEM, 
		Weighted_OSEM, Weighted_MLEM, Weighted_MBSREM, Weighted_BSREM, Weighted_ROSEM, Weighted_RBI, Weighted_COSEM, 
		TV_OSEM, TV_MLEM, TV_MBSREM, TV_BSREM, TV_ROSEM, TV_RBI, TV_COSEM, 
		AD_OSEM, AD_MLEM, AD_MBSREM, AD_BSREM, AD_ROSEM, AD_RBI, AD_COSEM, 
		APLS_OSEM, APLS_MLEM, APLS_MBSREM, APLS_BSREM, APLS_ROSEM, APLS_RBI, APLS_COSEM, 
		TGV_OSEM, TGV_MLEM, TGV_MBSREM, TGV_BSREM, TGV_ROSEM, TGV_RBI, TGV_COSEM,
		NLM_OSEM, NLM_MLEM, NLM_MBSREM, NLM_BSREM, NLM_ROSEM, NLM_RBI, NLM_COSEM,
		custom_OSEM, custom_MLEM, custom_MBSREM, custom_BSREM, custom_ROSEM, custom_RBI, custom_COSEM;
	af::array C_co = af::constant(0.f, 1, 1), C_aco = af::constant(0.f, 1, 1), C_osl = af::constant(0.f, 1, 1);
	af::array im_mlem, rhs_mlem, im_os, rhs_os, im_os_blurred, im_mlem_blurred;
} AF_im_vectors;

// Struct for the regularization parameters
typedef struct Beta_ {
	float MRP_OSEM, MRP_OSEMMAP, MRP_MLEM, MRP_MLEMMAP, MRP_MBSREM, MRP_BSREM, MRP_ROSEM, MRP_ROSEMOSL, MRP_RBI, MRP_RBIMAP, MRP_COSEM, MRP_COSEMMAP,
		Quad_OSEM, Quad_OSEMMAP, Quad_MLEM, Quad_MLEMMAP, Quad_MBSREM, Quad_BSREM, Quad_ROSEM, Quad_ROSEMOSL, Quad_RBI, Quad_RBIMAP, Quad_COSEM, Quad_COSEMMAP,
		Huber_OSEM, Huber_OSEMMAP, Huber_MLEM, Huber_MLEMMAP, Huber_MBSREM, Huber_BSREM, Huber_ROSEM, Huber_ROSEMOSL, Huber_RBI, Huber_RBIMAP, Huber_COSEM, Huber_COSEMMAP,
		L_OSEM, L_OSEMMAP, L_MLEM, L_MLEMMAP, L_MBSREM, L_BSREM, L_ROSEM, L_ROSEMOSL, L_RBI, L_RBIMAP, L_COSEM, L_COSEMMAP,
		FMH_OSEM, FMH_OSEMMAP, FMH_MLEM, FMH_MLEMMAP, FMH_MBSREM, FMH_BSREM, FMH_ROSEM, FMH_ROSEMOSL, FMH_RBI, FMH_RBIMAP, FMH_COSEM, FMH_COSEMMAP,
		Weighted_OSEM, Weighted_OSEMMAP, Weighted_MLEM, Weighted_MLEMMAP, Weighted_MBSREM, Weighted_BSREM, Weighted_ROSEM, Weighted_ROSEMOSL, Weighted_RBI, Weighted_RBIMAP, Weighted_COSEM, Weighted_COSEMMAP,
		TV_OSEM, TV_OSEMMAP, TV_MLEM, TV_MLEMMAP, TV_MBSREM, TV_BSREM, TV_ROSEM, TV_ROSEMOSL, TV_RBI, TV_RBIMAP, TV_COSEM, TV_COSEMMAP,
		AD_OSEM, AD_OSEMMAP, AD_MLEM, AD_MLEMMAP, AD_MBSREM, AD_BSREM, AD_ROSEM, AD_ROSEMOSL, AD_RBI, AD_RBIMAP, AD_COSEM, AD_COSEMMAP,
		APLS_OSEM, APLS_OSEMMAP, APLS_MLEM, APLS_MLEMMAP, APLS_MBSREM, APLS_BSREM, APLS_ROSEM, APLS_ROSEMOSL, APLS_RBI, APLS_RBIMAP, APLS_COSEM, APLS_COSEMMAP,
		TGV_OSEM, TGV_OSEMMAP, TGV_MLEM, TGV_MLEMMAP, TGV_MBSREM, TGV_BSREM, TGV_ROSEM, TGV_ROSEMOSL, TGV_RBI, TGV_RBIMAP, TGV_COSEM, TGV_COSEMMAP,
		NLM_OSEM, NLM_OSEMMAP, NLM_MLEM, NLM_MLEMMAP, NLM_MBSREM, NLM_BSREM, NLM_ROSEM, NLM_ROSEMOSL, NLM_RBI, NLM_RBIMAP, NLM_COSEM, NLM_COSEMMAP,
		custom_OSEM, custom_OSEMMAP, custom_MLEM, custom_MLEMMAP, custom_MBSREM, custom_BSREM, custom_ROSEM, custom_ROSEMOSL, custom_RBI, custom_RBIMAP, custom_COSEM, custom_COSEMMAP;
} Beta;

// Struct for various parameters, mainly various weights and coefficients
typedef struct Weighting_ {
	af::array tr_offsets, weights_quad, weights_TV, weights_huber, fmh_weights, a_L, weighted_weights, UU, Amin, D;
	af::array dU_OSEM, dU_MLEM, dU_BSREM, dU_MBSREM, dU_ROSEM, dU_RBI, dU_COSEM;
	af::array NLM_ref, gaussianNLM;
	float *lambda = nullptr, *lambda_MBSREM = nullptr, *lambda_BSREM = nullptr, *lambda_ROSEM = nullptr, *lambda_DRAMA = nullptr, h_ACOSEM = 1.f, TimeStepAD, KAD, w_sum = 0.f;
	float epsilon_mramla = 1e8f, U, NLM_gauss = 1.f, h2 = 1.f, huber_delta = 0.f, ACOSEM_rhs = 0.f, h_ACOSEM_2 = 1.f;
	uint32_t alku_fmh = 0u, mean_type = 0u;
	af_flux_function FluxType;
	af_diffusion_eq DiffusionType;
	uint32_t Ndx = 1u, Ndy = 1u, Ndz = 0u, NiterAD = 1u, dimmu, inffi, Nlx = 1u, Nly = 1u, Nlz = 0u;
	bool med_no_norm = false, MBSREM_prepass = false, NLM_MRP = false, NLTV = false, NLM_anatomical = false, deconvolution = false;
	uint32_t g_dim_x = 0u, g_dim_y = 0u, g_dim_z = 0u;
} Weighting;

// Struct for boolean operators indicating whether a certain method is selected
typedef struct RecMethods_ {
	bool MLEM = false, OSEM = false, MRAMLA = false, RAMLA = false, ROSEM = false, RBI = false, DRAMA = false, COSEM = false, ECOSEM = false, ACOSEM = false;
	bool MRP = false, Quad = false, Huber = false, L = false, FMH = false, WeightedMean = false, TV = false, AD = false, APLS = false, TGV = false, NLM = false;
	bool OSLMLEM = false, MAPMLEM = false, OSLOSEM = false, MAPOSEM = false, MBSREM = false, BSREM = false, ROSEMMAP = false, ROSEMOSL = false, RBIMAP = false, RBIOSL = false;
	bool MAP = false;
	bool CUSTOM = false;
	uint32_t OSLCOSEM = 0u, MAPCOSEM = 0u;
} RecMethods;

// MATLAB output arrays
typedef struct matlabArrays_ {
	mxArray* mlem, * osem, * ramla, * ramlaM, * rosem, * rbi, * drama, * cosem, * ecosem, * acosem,
		* mrp_mlem, * quad_mlem, * Huber_mlem, * L_mlem, * fmh_mlem, * weighted_mlem, * TV_mlem, * AD_mlem, * APLS_mlem, * TGV_mlem, * NLM_mlem,
		* mrp_osem, * quad_osem, * Huber_osem, * L_osem, * fmh_osem, * weighted_osem, * TV_osem, * AD_osem, * APLS_osem, * TGV_osem, * NLM_osem,
		* mrp_bsrem, * quad_bsrem, * Huber_bsrem, * L_bsrem, * fmh_bsrem, * weighted_bsrem, * TV_bsrem, * AD_bsrem, * APLS_bsrem, * TGV_bsrem, * NLM_bsrem,
		* mrp_mbsrem, * quad_mbsrem, * Huber_mbsrem, * L_mbsrem, * fmh_mbsrem, * weighted_mbsrem, * TV_mbsrem, * AD_mbsrem, * APLS_mbsrem, * TGV_mbsrem, * NLM_mbsrem,
		* mrp_rosem, * quad_rosem, * Huber_rosem, * L_rosem, * fmh_rosem, * weighted_rosem, * TV_rosem, * AD_rosem, * APLS_rosem, * TGV_rosem, * NLM_rosem,
		* mrp_rbi, * quad_rbi, * Huber_rbi, * L_rbi, * fmh_rbi, * weighted_rbi, * TV_rbi, * AD_rbi, * APLS_rbi, * TGV_rbi, * NLM_rbi,
		* mrp_cosem, * quad_cosem, * Huber_cosem, * L_cosem, * fmh_cosem, * weighted_cosem, * TV_cosem, * AD_cosem, * APLS_cosem, * TGV_cosem, * NLM_cosem,
		* custom_osem, * custom_mlem, * custom_bsrem, * custom_mbsrem, * custom_rosem, * custom_rbi, * custom_cosem;
	mxArray* c_osl_custom, *D_custom;
	float* ele_os, * ele_ml, * ele_ramla, * ele_ramlaM, * ele_rosem, * ele_rbi, * ele_drama, * ele_cosem, * ele_ecosem, * ele_acosem,
		* ele_mrp_mlem, * ele_quad_mlem, * ele_Huber_mlem, * ele_L_mlem, * ele_fmh_mlem, * ele_weighted_mlem, * ele_TV_mlem, * ele_AD_mlem, * ele_APLS_mlem, * ele_TGV_mlem, * ele_NLM_mlem,
		* ele_mrp_osem, * ele_quad_osem, * ele_Huber_osem, * ele_L_osem, * ele_fmh_osem, * ele_weighted_osem, * ele_TV_osem, * ele_AD_osem, * ele_APLS_osem, * ele_TGV_osem, * ele_NLM_osem,
		* ele_mrp_bsrem, * ele_quad_bsrem, * ele_Huber_bsrem, * ele_L_bsrem, * ele_fmh_bsrem, * ele_weighted_bsrem, * ele_TV_bsrem, * ele_AD_bsrem, * ele_TGV_bsrem, * ele_APLS_bsrem, * ele_NLM_bsrem,
		* ele_mrp_mbsrem, * ele_quad_mbsrem, * ele_Huber_mbsrem, * ele_L_mbsrem, * ele_fmh_mbsrem, * ele_weighted_mbsrem, * ele_TV_mbsrem, * ele_AD_mbsrem, * ele_TGV_mbsrem, * ele_APLS_mbsrem, * ele_NLM_mbsrem,
		* ele_mrp_rosem, * ele_quad_rosem, * ele_Huber_rosem, * ele_L_rosem, * ele_fmh_rosem, * ele_weighted_rosem, * ele_TV_rosem, * ele_AD_rosem, * ele_TGV_rosem, * ele_APLS_rosem, * ele_NLM_rosem,
		* ele_mrp_rbi, * ele_quad_rbi, * ele_Huber_rbi, * ele_L_rbi, * ele_fmh_rbi, * ele_weighted_rbi, * ele_TV_rbi, * ele_AD_rbi, * ele_TGV_rbi, * ele_APLS_rbi, * ele_NLM_rbi,
		* ele_mrp_cosem, * ele_quad_cosem, * ele_Huber_cosem, * ele_L_cosem, * ele_fmh_cosem, * ele_weighted_cosem, * ele_TV_cosem, * ele_AD_cosem, * ele_TGV_cosem, * ele_APLS_cosem, * ele_NLM_cosem,
		* ele_custom_osem, * ele_custom_mlem, * ele_custom_bsrem, * ele_custom_mbsrem, * ele_custom_rosem, * ele_custom_rbi, * ele_custom_cosem;
	float* ele_c_osl_custom, *ele_D_custom;
} matlabArrays;

typedef struct _kernelStruct {
#ifdef OPENCL
	cl::Kernel kernelNLM;
	cl::CommandQueue* af_queue;
#else
	CUfunction kernelNLM = NULL;
	CUstream* af_cuda_stream = nullptr;
#endif
} kernelStruct;

// Function for loading the data and forming the initial data variables (initial image estimates, etc.)
void form_data_variables(AF_im_vectors &vec, Beta &beta, Weighting &w_vec, const mxArray* options, const uint32_t Nx, const uint32_t Ny,
	const uint32_t Nz, const uint32_t Niter, const af::array &x0, const uint32_t im_dim, const size_t koko_l, const RecMethods &MethodList, TVdata &data, 
	const uint32_t subsets, const uint32_t osa_iter0, const bool use_psf, const bool saveIter, const uint32_t Nt);

// Get the reconstruction methods used
void get_rec_methods(const mxArray *options, RecMethods &MethodList);

// MATLAB output
void create_matlab_output(matlabArrays &ArrayList, const mwSize *dimmi, const RecMethods &MethodList, const uint32_t dim_n);

// Transfer device data back to host MATLAB cell
void device_to_host_cell(matlabArrays &ArrayList, const RecMethods &MethodList, AF_im_vectors & vec, uint32_t &oo, mxArray *cell, Weighting& w_vec);

// Compute the epsilon value for the MBSREM/MRAMLA
float MBSREM_epsilon(const af::array &Sino, const float epps, const uint32_t randoms_correction, const af::array& randoms, const af::array& D, 
	const bool TOF, const int64_t nBins);

// Batch functions (scalar and vector or vector and matrix)
af::array batchMinus(const af::array &lhs, const af::array &rhs);

af::array batchPlus(const af::array &lhs, const af::array &rhs);

af::array batchMul(const af::array &lhs, const af::array &rhs);

af::array batchDiv(const af::array &lhs, const af::array &rhs);

af::array batchNotEqual(const af::array &lhs, const af::array &rhs);

// Create a zero-padded image
af::array padding(const af::array& im, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, 
	const bool zero_pad = false);

// Reconstruction methods
//af::array MLEM(const af::array &im, const af::array &Summ, const af::array &rhs);

//af::array OSL_MLEM(const af::array &im, const af::array &Summ, const af::array &rhs, const af::array &dU, const float beta);

af::array EM(const af::array &im, const af::array &Summ, const af::array &rhs);

af::array OSL(const af::array& Summ, const af::array& dU, const float beta, const float epps);

af::array MBSREM(const af::array & im, const af::array & rhs, const float U, const af::array & pj3, const float* lam, const uint32_t iter, const uint32_t im_dim,
	const float beta, const af::array &dU, const af::array & Summ, const float epps);

af::array BSREM(const af::array &im, const af::array &rhs, const float *lam, const uint32_t iter);

//af::array COSEM(const af::array &im, const af::array &C_co, const af::array &D);

af::array ECOSEM(const af::array &im, const af::array &D, const af::array &OSEM_apu, const af::array &COSEM_apu, const float epps);

//af::array ACOSEM(const af::array & im, const af::array & C_aco, const af::array & D, const float h);

af::array ROSEM(const af::array &im, const af::array &Summ, const af::array &rhs, const float *lam, const uint32_t iter);

af::array RBI(const af::array& im, const af::array& Summ, const af::array& rhs, const af::array& D = af::constant(0.f, 1, 1), const float beta = 0.f, const af::array& dU = af::constant(0.f, 1, 1));

af::array DRAMA(const af::array &im, const af::array &Summ, const af::array &rhs, const float *lam, const uint32_t iter, const uint32_t sub_iter,
	const uint32_t subsets);

af::array MAP(const af::array &im, const float lam, const float beta, const af::array &dU, const float epps);

af::array COSEM(const af::array &im, const af::array &C_co, const af::array &D, const float h, const uint32_t COSEM_TYPE);

// Priors
af::array MRP(const af::array &im, const uint32_t medx, const uint32_t medy, const uint32_t medz, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, 
	const float epps, const af::array &offsets, const bool med_no_norm, const uint32_t im_dim);

af::array Quadratic_prior(const af::array& im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const uint32_t Nx, const uint32_t Ny,
	const uint32_t Nz, const uint32_t inffi, const af::array& offsets, const af::array& weights_quad, const uint32_t im_dim);

af::array Huber_prior(const af::array& im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const uint32_t Nx, const uint32_t Ny,
	const uint32_t Nz, const uint32_t inffi, const af::array& offsets, const af::array& weights_huber, const uint32_t im_dim, const float delta);

af::array FMH(const af::array &im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, 
	const float epps, const uint32_t inffi, const af::array &offsets, const af::array &fmh_weights, const bool med_no_norm, const uint32_t alku_fmh, 
	const uint32_t im_dim);

af::array L_filter(const af::array &im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, 
	const float epps, const af::array &offsets, const af::array &a_L, const bool med_no_norm, const uint32_t im_dim);

af::array Weighted_mean(const af::array &im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const uint32_t Nx, const uint32_t Ny, 
	const uint32_t Nz, const float epps, const af::array &weighted_weights, const bool med_no_norm, const uint32_t im_dim, 
	const uint32_t mean_type, const float w_sum);

af::array AD(const af::array &im, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const float epps, const float TimeStepAD, const float KAD, 
	const uint32_t NiterAD, const af_flux_function FluxType, const af_diffusion_eq DiffusionType, const bool med_no_norm);

af::array TVprior(const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const TVdata &S, const af::array& im, const float epps, const uint32_t TVtype, 
	const Weighting & w_vec, const af::array& offsets);

af::array TGV(const af::array &im, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const uint32_t maxits, const float alpha, const float beta);

//af::array NLM(const af::array& im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const uint32_t Nlx, const uint32_t Nly, const uint32_t Nlz,
//	const float h2, const float epps, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const bool NLM_anatomical, const af::array& gaussianNLM,
//	const bool NLTV, const bool NLM_MRP, const af::array& NLM_ref);

void reconstruction_AF_matrixfree(const size_t koko, const uint16_t* lor1, const float* z_det, const float* x, const float* y, const mxArray* Sin,
	const mxArray* sc_ra, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const uint32_t Niter, const mxArray* options, const float dx,
	const float dy, const float dz, const float bx, const float by, const float bz, const float bzb, const float maxxx, const float maxyy, const float zmax,
	const float NSlices, const int64_t* pituus, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t size_x, const uint32_t TotSinos,
	mxArray* cell, const mwSize* dimmi, const bool verbose, const uint32_t randoms_correction, const uint32_t attenuation_correction,
	const uint32_t normalization, const float* atten, const size_t size_atten, const float* norm, const size_t size_norm, const uint32_t subsets,
	const float epps, const char* k_path, const uint32_t Nt, const uint32_t* pseudos, const uint32_t det_per_ring, const uint32_t prows, const uint16_t* L,
	const uint8_t raw, const size_t size_z, const bool osem_bool, const char* fileName, const bool use_psf, const float tube_width,
	const float crystal_size_z, const float* x_center, const float* y_center, const float* z_center, const size_t size_center_x, const size_t size_center_y,
	const size_t size_of_x, const size_t size_center_z, const uint32_t projector_type, const char* header_directory, const bool precompute,
	const uint32_t device, const int32_t dec, const uint16_t n_rays, const uint16_t n_rays3D, const float cr_pz, const bool use_64bit_atomics, uint32_t n_rekos,
	const uint32_t n_rekos_mlem, const uint8_t* reko_type, const uint8_t* reko_type_mlem, const float global_factor, const float bmin, const float bmax,
	const float Vmax, const float* V, const size_t size_V, const float* gaussian, const size_t size_gauss, const bool saveIter, const bool TOF, const int64_t TOFSize,
	const float sigma_x, const float* TOFCenter, const int64_t nBins);

//void reconstruction_AF_matrixfree_CUDA(const size_t koko, const uint16_t* lor1, const float* z_det, const float* x, const float* y, const mxArray* Sin,
//	const mxArray* sc_ra, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const uint32_t Niter, const mxArray* options, const float dx,
//	const float dy, const float dz, const float bx, const float by, const float bz, const float bzb, const float maxxx, const float maxyy, const float zmax,
//	const float NSlices, const uint32_t* pituus, const uint32_t* xy_index, const uint16_t* z_index, const uint32_t size_x, const uint32_t TotSinos,
//	mxArray* cell, const mwSize* dimmi, const bool verbose, const uint32_t randoms_correction, const uint32_t attenuation_correction,
//	const uint32_t normalization, const float* atten, const size_t size_atten, const float* norm, const size_t size_norm, const uint32_t subsets,
//	const float epps, const char* k_path, const uint32_t Nt, const uint32_t* pseudos, const uint32_t det_per_ring, const uint32_t prows, const uint16_t* L,
//	const uint8_t raw, const size_t size_z, const bool osem_bool, const char* fileName, const bool force_build, const float tube_width,
//	const float crystal_size_z, const float* x_center, const float* y_center, const float* z_center, const size_t size_center_x, const size_t size_center_y,
//	const size_t size_of_x, const size_t size_center_z, const uint32_t projector_type, const char* header_directory, const bool precompute,
//	const uint32_t device, const int32_t dec, const uint16_t n_rays, const uint16_t n_rays3D, const float cr_pz, const bool use_64bit_atomics, uint32_t n_rekos,
//	const uint32_t n_rekos_mlem, const uint8_t* reko_type, const uint8_t* reko_type_mlem, const float global_factor, const float bmin, const float bmax,
//	const float Vmax, const float* V, const size_t size_V, const float* gaussian, const size_t size_gauss);

//af::array im2col_3D(const af::array& A, const uint32_t blocksize1, const uint32_t blocksize2, const uint32_t blocksize3);
//
//af::array sparseSum(const af::array& W, const uint32_t s);

af::array computeConvolution(const af::array& vec, const af::array& g, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const Weighting& w_vec,
	const uint32_t n_rekos);

void deblur(af::array& vec, const af::array& g, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const Weighting& w_vec, const uint32_t iter, 
	const uint32_t subsets, const float epps, const bool saveIter);

void computeDeblur(AF_im_vectors& vec, const af::array& g, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const Weighting& w_vec,
	const RecMethods& MethodList, const uint32_t iter, const uint32_t subsets, const float epps, const bool saveIter);

void computeDeblurMLEM(AF_im_vectors& vec, const af::array& g, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const Weighting& w_vec,
	const RecMethods& MethodList, const uint32_t iter, const uint32_t subsets, const float epps, const bool saveIter);

void computeOSEstimatesIter(AF_im_vectors& vec, Weighting& w_vec, const RecMethods& MethodList, const uint32_t im_dim, const float epps,
	const uint32_t iter, const uint32_t osa_iter0, const uint32_t subsets, const Beta& beta, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz,
	const TVdata& data, const uint32_t n_rekos2, const kernelStruct& OpenCLStruct, const bool saveIter);

void computeMLEstimates(AF_im_vectors& vec, Weighting& w_vec, const RecMethods& MethodList, const uint32_t im_dim, const float epps,
	const uint32_t iter, const uint32_t subsets, const Beta& beta, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz,
	const TVdata& data, const af::array& Summ_mlem, bool& break_iter, const kernelStruct& OpenCLStruct, const bool saveIter);

af::array NLM(const af::array& im, Weighting& w_vec, const float epps, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const kernelStruct& OpenCLStruct);
