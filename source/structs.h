#pragma once
#define DEBUG true
#include <arrayfire.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include "precomp.h"
#ifdef OPENCL
#include <af/opencl.h>
#else
#include <af/cuda.h>
#endif


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
	af::array C_co = af::constant(0.f, 1, 1), C_aco = af::constant(0.f, 1, 1), C_osl = af::constant(0.f, 1, 1);
	af::array wLSQR = af::constant(0.f, 1, 1), fLSQR = af::constant(0.f, 1, 1);
	af::array im_os, rhs_os, im_os_blurred, meanFP, meanBP;
} AF_im_vectors;

// Struct for various parameters, mainly various weights and coefficients
typedef struct Weighting_ {
	af::array tr_offsets, weights_quad, weights_TV, weights_huber, fmh_weights, a_L, weighted_weights, UU, Amin, D, weights_RDP;
	std::vector<af::array> dU;
	af::array gaussianNLM, gFilter;
	float* lambda = nullptr, * lambda_MBSREM = nullptr, * lambda_BSREM = nullptr, * lambda_ROSEM = nullptr, * lambda_DRAMA = nullptr, h_ACOSEM = 1.f, TimeStepAD, KAD, w_sum = 0.f,
		* lambda_PKMA = nullptr, * alpha_PKMA = nullptr, * sigma_PKMA = nullptr, * uv = nullptr, * angles = nullptr, * NLM_ref = nullptr;
	uint32_t* rekot = nullptr, * distInt = nullptr;
	uint8_t* maskFP = nullptr, * maskBP = nullptr;
	float epsilon_mramla = 1e8f, U, NLM_gauss = 1.f, h2 = 1.f, huber_delta = 0.f, ACOSEM_rhs = 0.f, h_ACOSEM_2 = 1.f, RDP_gamma, tauCP = 0.f,
		sigmaCP = 0.f, thetaCP = 0.f, kerroin4 = 0.f, sigma2CP = 0.f, dPitchX, dPitchY, betaLSQR = 0.f, alphaLSQR = 0.f, thetaLSQR = 0.f, rhoLSQR = 0.f, 
		phiLSQR = 0.f;
	uint32_t alku_fmh = 0u, mean_type = 0u, deblur_iterations = 0;
	af_flux_function FluxType;
	af_diffusion_eq DiffusionType;
	uint32_t Ndx = 1u, Ndy = 1u, Ndz = 0u, NiterAD = 1u, dimmu, inffi, Nlx = 1u, Nly = 1u, Nlz = 0u;
	bool med_no_norm = false, MBSREM_prepass = false, NLM_MRP = false, NLTV = false, NLM_anatomical = false, deconvolution = false;
	uint32_t g_dim_x = 0u, g_dim_y = 0u, g_dim_z = 0u;
	uint32_t size_y = 0U, size_x = 0U;
	int64_t nProjections = 0LL;
	//cl_float2 dPitch = { 0.f, 0.f };
	uint32_t nPriors = 0U, nMAP = 0U, nMAPML = 0U, nMLEM = 0U, nOS = 0U, nTot = 0U, nMAPOS = 0U, nPriorsTot = 0U;
	std::vector<int32_t> mIt;
} Weighting;

// Struct for boolean operators indicating whether a certain method is selected
typedef struct RecMethods_ {
	bool OSEM = false, MRAMLA = false, RAMLA = false, ROSEM = false, RBI = false, DRAMA = false, COSEM = false, ECOSEM = false,
		ACOSEM = false, LSQR = false, CGLS = false, CPLS = false;
	bool MRP = false, Quad = false, Huber = false, L = false, FMH = false, WeightedMean = false, TV = false, AD = false, APLS = false, TGV = false, NLM = false, RDP = false;
	bool OSLOSEM = false, MAPOSEM = false, MBSREM = false, BSREM = false, ROSEMMAP = false, ROSEMOSL = false, RBIMAP = false, RBIOSL = false,
		PKMA = false, CPTV = false;
	bool MAP = false;
	bool CUSTOM = false;
	uint32_t OSLCOSEM = 0u, MAPCOSEM = 0u;
} RecMethods;

#ifdef OPENCL
// Struct for the various estimates used in OpenCL kernels
typedef struct _OpenCL_im_vectors {
	cl::Buffer d_rhs_os, d_meanFP, d_meanBP;
	//cl::Image3D d_image_os, d_image_os_int;
	cl::Image3D d_image_os, d_image_os_int;
} OpenCL_im_vectors;
#else

#endif