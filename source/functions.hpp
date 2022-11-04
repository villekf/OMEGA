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
#include "ProjectorClass.h"
//#include "mexFunktio.h"
#ifdef OPENCL
#else
#include <nvrtc.h>
#include <cublas.h>
#include <cuda.h>
#include <string>
#include <iostream>
#include <fstream>

const char* getErrorString(CUresult error);

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(CUresult code, const char* file, int line, bool abort = true)
{
	if (code != CUDA_SUCCESS)
	{
		const char* errstr;
		cuGetErrorString(code, &errstr);
		mexPrintf("GPUassert: %s %s %d\n", errstr, file, line);
	}
}
#endif
//#include "ProjectorClass.h"

#pragma pack(1) 
#pragma warning(disable : 4996)

// MATLAB output arrays
//typedef struct matlabArrays_ {
//	mxArray* mlem, * osem, * ramla, * ramlaM, * rosem, * rbi, * drama, * cosem, * ecosem, * acosem,
//		* mrp_mlem, * quad_mlem, * Huber_mlem, * L_mlem, * fmh_mlem, * weighted_mlem, * TV_mlem, * AD_mlem, * APLS_mlem, * TGV_mlem, * NLM_mlem,
//		* mrp_osem, * quad_osem, * Huber_osem, * L_osem, * fmh_osem, * weighted_osem, * TV_osem, * AD_osem, * APLS_osem, * TGV_osem, * NLM_osem,
//		* mrp_bsrem, * quad_bsrem, * Huber_bsrem, * L_bsrem, * fmh_bsrem, * weighted_bsrem, * TV_bsrem, * AD_bsrem, * APLS_bsrem, * TGV_bsrem, * NLM_bsrem,
//		* mrp_mbsrem, * quad_mbsrem, * Huber_mbsrem, * L_mbsrem, * fmh_mbsrem, * weighted_mbsrem, * TV_mbsrem, * AD_mbsrem, * APLS_mbsrem, * TGV_mbsrem, * NLM_mbsrem,
//		* mrp_rosem, * quad_rosem, * Huber_rosem, * L_rosem, * fmh_rosem, * weighted_rosem, * TV_rosem, * AD_rosem, * APLS_rosem, * TGV_rosem, * NLM_rosem,
//		* mrp_rbi, * quad_rbi, * Huber_rbi, * L_rbi, * fmh_rbi, * weighted_rbi, * TV_rbi, * AD_rbi, * APLS_rbi, * TGV_rbi, * NLM_rbi,
//		* mrp_cosem, * quad_cosem, * Huber_cosem, * L_cosem, * fmh_cosem, * weighted_cosem, * TV_cosem, * AD_cosem, * APLS_cosem, * TGV_cosem, * NLM_cosem,
//		* custom_osem, * custom_mlem, * custom_bsrem, * custom_mbsrem, * custom_rosem, * custom_rbi, * custom_cosem;
//	mxArray* c_osl_custom, *D_custom;
//	float* ele_os, * ele_ml, * ele_ramla, * ele_ramlaM, * ele_rosem, * ele_rbi, * ele_drama, * ele_cosem, * ele_ecosem, * ele_acosem,
//		* ele_mrp_mlem, * ele_quad_mlem, * ele_Huber_mlem, * ele_L_mlem, * ele_fmh_mlem, * ele_weighted_mlem, * ele_TV_mlem, * ele_AD_mlem, * ele_APLS_mlem, * ele_TGV_mlem, * ele_NLM_mlem,
//		* ele_mrp_osem, * ele_quad_osem, * ele_Huber_osem, * ele_L_osem, * ele_fmh_osem, * ele_weighted_osem, * ele_TV_osem, * ele_AD_osem, * ele_APLS_osem, * ele_TGV_osem, * ele_NLM_osem,
//		* ele_mrp_bsrem, * ele_quad_bsrem, * ele_Huber_bsrem, * ele_L_bsrem, * ele_fmh_bsrem, * ele_weighted_bsrem, * ele_TV_bsrem, * ele_AD_bsrem, * ele_TGV_bsrem, * ele_APLS_bsrem, * ele_NLM_bsrem,
//		* ele_mrp_mbsrem, * ele_quad_mbsrem, * ele_Huber_mbsrem, * ele_L_mbsrem, * ele_fmh_mbsrem, * ele_weighted_mbsrem, * ele_TV_mbsrem, * ele_AD_mbsrem, * ele_TGV_mbsrem, * ele_APLS_mbsrem, * ele_NLM_mbsrem,
//		* ele_mrp_rosem, * ele_quad_rosem, * ele_Huber_rosem, * ele_L_rosem, * ele_fmh_rosem, * ele_weighted_rosem, * ele_TV_rosem, * ele_AD_rosem, * ele_TGV_rosem, * ele_APLS_rosem, * ele_NLM_rosem,
//		* ele_mrp_rbi, * ele_quad_rbi, * ele_Huber_rbi, * ele_L_rbi, * ele_fmh_rbi, * ele_weighted_rbi, * ele_TV_rbi, * ele_AD_rbi, * ele_TGV_rbi, * ele_APLS_rbi, * ele_NLM_rbi,
//		* ele_mrp_cosem, * ele_quad_cosem, * ele_Huber_cosem, * ele_L_cosem, * ele_fmh_cosem, * ele_weighted_cosem, * ele_TV_cosem, * ele_AD_cosem, * ele_TGV_cosem, * ele_APLS_cosem, * ele_NLM_cosem,
//		* ele_custom_osem, * ele_custom_mlem, * ele_custom_bsrem, * ele_custom_mbsrem, * ele_custom_rosem, * ele_custom_rbi, * ele_custom_cosem;
//	float* ele_c_osl_custom, *ele_D_custom;
//} matlabArrays;

//typedef struct _kernelStruct {
//#ifdef OPENCL
//	cl::Kernel kernelNLM;
//	cl::Kernel kernelMed;
//	cl::CommandQueue* af_queue;
//#else
//	CUfunction kernelNLM = NULL;
//	CUfunction kernelMed = NULL;
//	CUstream* af_cuda_stream = nullptr;
//#endif
//} kernelStruct;

// Function for loading the data and forming the initial data variables (initial image estimates, etc.)
void form_data_variables(AF_im_vectors& vec, std::vector<float>& beta, Weighting& w_vec, const mxArray* options, scalarStruct& inputScalars,
	const RecMethods& MethodList, TVdata& data, const uint32_t Nt, const uint32_t iter0, std::vector<std::vector<float*>>& imEstimates);

// Get the reconstruction methods used
void get_rec_methods(const mxArray *options, RecMethods &MethodList);

// MATLAB output
//void create_matlab_output(matlabArrays &ArrayList, const mwSize *dimmi, const RecMethods &MethodList, const uint32_t dim_n);

// Transfer device data back to host MATLAB cell
void device_to_host_cell(const RecMethods &MethodList, AF_im_vectors & vec, uint32_t &oo, mxArray *cell, Weighting& w_vec,
	const mwSize* dimmi, const uint32_t dim_n, const std::vector<std::vector<float*>>& imEstimates, const scalarStruct& inputScalars);

// Compute the epsilon value for the MBSREM/MRAMLA
float MBSREM_epsilon(const af::array &Sino, const float epps, const uint32_t randoms_correction, const af::array& randoms, const af::array& D, 
	const bool TOF, const int64_t nBins, const bool CT = false);

// Batch functions (scalar and vector or vector and matrix)
af::array batchMinus(const af::array &lhs, const af::array &rhs);

af::array batchPlus(const af::array &lhs, const af::array &rhs);

af::array batchMul(const af::array &lhs, const af::array &rhs);

af::array batchDiv(const af::array &lhs, const af::array &rhs);

af::array batchNotEqual(const af::array &lhs, const af::array &rhs);

// Create a zero-padded image
af::array padding(const af::array& im, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz,
	const bool zero_pad = false, const uint32_t Nw = 1);

// Reconstruction methods
//af::array MLEM(const af::array &im, const af::array &Summ, const af::array &rhs);

//af::array OSL_MLEM(const af::array &im, const af::array &Summ, const af::array &rhs, const af::array &dU, const float beta);

af::array EM(const af::array &im, const af::array &Summ, const af::array &rhs);

af::array OSL(const af::array& Summ, const af::array& dU, const float beta, const float epps);

af::array MBSREM(const af::array& im, const af::array& rhs, const float U, const af::array& D, const float* lam, const uint32_t iter,
	const float beta, const af::array& dU, const af::array& Summ, const scalarStruct inputScalars);

af::array BSREM(const af::array &im, const af::array &rhs, const float *lam, const uint32_t iter, const af::array& Summ);

//af::array COSEM(const af::array &im, const af::array &C_co, const af::array &D);

af::array ECOSEM(const af::array &im, const af::array &D, const af::array &OSEM_apu, const af::array &COSEM_apu, const float epps);

//af::array ACOSEM(const af::array & im, const af::array & C_aco, const af::array & D, const float h);

af::array ROSEM(const af::array &im, const af::array &Summ, const af::array &rhs, const float *lam, const uint32_t iter);

af::array RBI(const af::array& im, const af::array& Summ, const af::array& rhs, const af::array& D = af::constant(0.f, 1, 1), const float beta = 0.f, const af::array& dU = af::constant(0.f, 1, 1));

af::array DRAMA(const af::array &im, const af::array &Summ, const af::array &rhs, const float *lam, const uint32_t iter, const uint32_t sub_iter,
	const uint32_t subsets);

af::array MAP(const af::array &im, const float lam, const float beta, const af::array &dU, const float epps);

af::array COSEM(const af::array &im, const af::array &C_co, const af::array &D, const float h, const uint32_t COSEM_TYPE);

af::array PKMA(const af::array& im, const af::array& Summ, const af::array& rhs, const float* lam, const float* alpha, const float* sigma, const af::array& D,
	const uint32_t iter, const uint32_t osa_iter, const uint32_t subsets, const float epps, const float beta, const af::array& dU);

void LSQR(af::array& im, const af::array& rhs, const scalarStruct& inputScalars, Weighting& w_vec, const uint32_t iter, AF_im_vectors& vec);

// Priors
af::array MRP(const af::array &im, const uint32_t medx, const uint32_t medy, const uint32_t medz, const scalarStruct& inputScalars, 
	const af::array &offsets, const bool med_no_norm, ProjectorClass& proj);

af::array Quadratic_prior(const af::array& im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const scalarStruct& inputScalars, 
	const uint32_t inffi, const af::array& offsets, const af::array& weights_quad);

af::array Huber_prior(const af::array& im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const scalarStruct& inputScalars, 
	const uint32_t inffi, const af::array& offsets, const af::array& weights_huber, const float delta);

af::array FMH(const af::array &im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const scalarStruct& inputScalars,
	const uint32_t inffi, const af::array &offsets, const af::array &fmh_weights, const bool med_no_norm, const uint32_t alku_fmh);

af::array L_filter(const af::array &im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const scalarStruct& inputScalars, 
	const af::array &offsets, const af::array &a_L, const bool med_no_norm);

af::array Weighted_mean(const af::array &im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const scalarStruct& inputScalars, 
	const af::array &weighted_weights, const bool med_no_norm, const uint32_t mean_type, const float w_sum);

af::array AD(const af::array &im, const scalarStruct& inputScalars, const float TimeStepAD, const float KAD,
	const uint32_t NiterAD, const af_flux_function FluxType, const af_diffusion_eq DiffusionType, const bool med_no_norm);

af::array TVprior(const scalarStruct& inputScalars, const TVdata &S, const af::array& im, const uint32_t TVtype, const Weighting & w_vec, 
	const af::array& offsets);

af::array TGV(const af::array &im, const scalarStruct& inputScalars, const uint32_t maxits, const float alpha, const float beta);

af::array RDP(const af::array& im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const scalarStruct& inputScalars, const af::array& weights_RDP,
	const float gamma, const af::array& offsets, const uint32_t inffi);

//af::array NLM(const af::array& im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const uint32_t Nlx, const uint32_t Nly, const uint32_t Nlz,
//	const float h2, const float epps, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const bool NLM_anatomical, const af::array& gaussianNLM,
//	const bool NLTV, const bool NLM_MRP, const af::array& NLM_ref);

void reconstruction_AF_matrixfree(const uint16_t* lor1, const float* z_det, const float* x, const float* y, const mxArray* Sin,
	const mxArray* sc_ra, scalarStruct inputScalars, const mxArray* options, const int64_t* pituus, const uint32_t* xy_index,
	const uint16_t* z_index, mxArray* cell, const mwSize* dimmi, const float* atten, const float* norm, const char* k_path,
	const uint32_t Nt, const uint32_t* pseudos, const uint32_t prows, const uint16_t* L, const char* fileName, const float* x_center,
	const float* y_center, const float* z_center, const char* header_directory, const uint32_t device, uint32_t n_rekos,
	const uint8_t* reko_type, const float* V, const float* gaussian, const size_t size_gauss, const float* TOFCenter);

af::array computeConvolution(const af::array& vec, const af::array& g, scalarStruct& inputScalars, const Weighting& w_vec,
	const uint32_t nRekos = 1);

void deblur(af::array& vec, const af::array& g, const scalarStruct& inputScalars, const Weighting& w_vec);

//void computeDeblur(af::array& vec, const af::array& g, const scalarStruct& inputScalars, const Weighting& w_vec,
//	const RecMethods& MethodList);

//void computeDeblurMLEM(AF_im_vectors& vec, const af::array& g, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const Weighting& w_vec,
//	const RecMethods& MethodList, const uint32_t iter, const uint32_t subsets, const float epps, const bool saveIter);

void computeOSEstimatesIter(AF_im_vectors& vec, Weighting& w_vec, const RecMethods& MethodList, const scalarStruct& inputScalars, const uint32_t iter,
	const std::vector<float>& beta, const TVdata& data, std::vector<std::vector<float*>>& imEstimates, ProjectorClass& proj, const af::array& g);

//void computeMLEstimates(AF_im_vectors& vec, Weighting& w_vec, const RecMethods& MethodList, const uint32_t im_dim, const float epps,
//	const uint32_t iter, const uint32_t subsets, const std::vector<float>& beta, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz,
//	const TVdata& data, const af::array& Summ_mlem, bool& break_iter, const kernelStruct& OpenCLStruct, const bool saveIter);

af::array NLM(ProjectorClass& proj, const af::array& im, Weighting& w_vec, const scalarStruct& inputScalars);

void computeForwardStep(const RecMethods& MethodList, af::array& y, af::array& input, const int64_t length, const scalarStruct& inputScalars,
	Weighting& w_vec, const af::array& randomsData);

void computeIntegralImage(const scalarStruct& inputScalars, const Weighting& w_vec, const int64_t length, af::array& outputFP, af::array& meanBP);

void computeOSEstimates(AF_im_vectors& vec, Weighting& w_vec, const RecMethods& MethodList, af::array* testi, const uint32_t iter,
	const uint32_t osa_iter, const scalarStruct& inputScalars, const std::vector<float>& beta, const TVdata& data, std::vector<int64_t>& length,
	bool& break_iter, const int64_t* pituus, std::vector<af::array>& Summ, af::array& E, const af::array& g, const af::array& D, const mxArray* Sin,
	ProjectorClass& proj);

void initializeRHS(AF_im_vectors& vec, scalarStruct& inputScalars);

void initializationStep(Weighting& w_vec, af::array& mData, AF_im_vectors& vec, ProjectorClass& proj, scalarStruct& inputScalars,
	std::vector<int64_t> length, uint64_t m_size, uint64_t st, const RecMethods& MethodList, uint32_t curIter, const af::array& g, 
	af::array& meanBP, std::vector<af::array>& Summ);

void forwardProjectionSPECT(af::array& fProj, const Weighting& w_vec, AF_im_vectors& vec, const scalarStruct& inputScalars,
	const int64_t length, const uint32_t uu);

void backprojectionSPECT(af::array& fProj, std::vector<af::array>& Summ, const Weighting& w_vec, AF_im_vectors& vec,
	const scalarStruct& inputScalars, const int64_t length, const uint32_t uu, const uint32_t osa_iter, const uint32_t iter,
	const uint8_t compute_norm_matrix, const uint32_t iter0);