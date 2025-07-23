#pragma once
// If true, displays a lot of debug messages
#define DEBUG false
// If true, uses the (correct) method to compute the inverse of measurement-based filtered PDHG, if false uses the same version as non-filtered one
// Only affects measurement-based filtering
#define FINVERSE true
// If false, uses the original EM preconditioner. If true, uses an EM preconditioner with the measurement data instead of ones.
#define SENSSCALE false
// If true, outputs all the multi-resolution volumes (and the reconstruction itself) as different cell elements
#define CELL false
#include <algorithm>
#include <vector>
#include <cmath>
#include <numeric>
#include "precomp.h"
#ifdef AF
#include <arrayfire.h>
#ifdef CUDA
#define AF_DEFINE_CUDA_TYPES 1
#include <af/cuda.h>
#include <cstring>
#elif defined(OPENCL)
#include <af/opencl.h>
#endif
#endif

// Used in ProjectorClass.h and ProjectorClassCUDA.h
#if defined(OPENCL)
#define OCL_CHECK(STATUS, MSG, RETURN) do { \
    if ((STATUS) != CL_SUCCESS) { \
        getErrorString(STATUS); \
        mexPrint(MSG); \
        return RETURN; \
    } \
} while(0)
#elif defined(CUDA)
#define CUDA_CHECK(STATUS, MSG, RETURN) do { \
    if ((STATUS) != CUDA_SUCCESS) { \
        getErrorString(STATUS); \
        mexPrint(MSG); \
        return RETURN; \
    } \
} while(0)
#endif

// Struct for the TV-prior
typedef struct TVdata_ {
#ifdef AF
	af::array refIm;
#endif
	bool TV_use_anatomical = false;
	float tau, TVsmoothing, T, C, eta, APLSsmoothing, TGVAlpha, TGVBeta, SATVPhi = 0.f;
	uint32_t TVtype = 0;
	uint32_t NiterTGV;
} TVdata;

#ifdef AF
// Struct for the various ArrayFire arrays
typedef struct AF_im_vectors_ {
	af::array C_co, dU;
	af::array rCGLS, meanFP, meanBP, apu, pPrevCP, p0CP2, fpCP2, p0CP, adapTypeA;
	std::vector<af::array> im_os, im_os_blurred, rhs_os, pCP, qProxTGV, vProxTGV, qProxTV, qProx, SAGASum;
	std::vector<af::array> wLSQR, fLSQR, uCP, uFISTA, fCGLS, rhsCP, fpCP, f0POCS, gradBB, imBB;
	std::vector<std::vector<af::array>> Summ, stochasticHelper;
} AF_im_vectors;
#endif

// Struct for various parameters, mainly various weights and coefficients
typedef struct Weighting_ {
#ifdef AF
	af::array tr_offsets, weights_quad, weights_TV, weights_huber, fmh_weights, a_L, weighted_weights, UU, Amin, weights_RDP, RDPref;
	std::vector<af::array> dU, M, preRef, dP, D, gradF;
	af::array gaussianNLM, gFilter, filter, filterIm, Ffilter;
	af_flux_function FluxType;
	af_diffusion_eq DiffusionType;
#endif
	TVdata data;
	float* lambda = nullptr, * alphaM = nullptr, * sigma_PKMA = nullptr, * uv = nullptr, * angles = nullptr, * NLM_ref = nullptr, * refImage = nullptr,
		* alphaPrecond = nullptr, * thetaCP = nullptr, * tauCP2 = nullptr, * tauCP = nullptr, * sigmaCP = nullptr, * sigma2CP = nullptr, * kerroin4 = nullptr, * lambdaFiltered = nullptr, 
		* weights = nullptr, *listCoord = nullptr, * RDP_ref = nullptr;
	uint32_t* rekot = nullptr;
    int32_t *distInt = nullptr, *distInt2 = nullptr;
	uint8_t* maskFP = nullptr, * maskBP = nullptr, * eFOVIndices = nullptr, *maskPrior = nullptr, *maskOffset = nullptr, *TOFIndices = nullptr;
	uint16_t* axIndex = nullptr, * trIndex = nullptr;
	float epsilon_mramla = 0.f, U = 1000000.f, h_ACOSEM = 1.f, TimeStepAD, KAD, w_sum = 0.f, h2 = 1.f, huber_delta = 0.f, ACOSEM_rhs = 0.f, h_ACOSEM_2 = 1.f, RDP_gamma = 1.f,
		dPitchX, dPitchY, betaLSQR = 0.f, alphaLSQR = 0.f, thetaLSQR = 0.f, rhoLSQR = 0.f, phiLSQR = 0.f, gammaCGLS = 0.f, gammaTempCGLS = 0.f, alphaCGLS = 0.f, nuIEM = 0.f, alphaCPTV = 1.f,
		gradV1 = 0.f, gradV2 = 0.f, betaReg = 0.f, alpha0CPTGV = 1.f, alpha1CPTGV = 1.f, betaFISTA = 1.f, tFISTA = 1.f, tNFista = 1.f, GGMRF_p = 0.f, GGMRF_q = 0.f, GGMRF_c = 0.f, GGMRF_pqc = 0.f,
		beta = 0.f, dtvg = 0.f, alphaPOCS = 0.2f, rMaxPOCS = 0.95f, POCSepps = 1e-4f, POCSalphaRed = 0.95f, NLAdaptiveConstant = 1e-5f;
	uint32_t alku_fmh = 0u, mean_type = 0u, powerIterations = 0, derivType = 0, gradInitIter = 0, filterIter = 0, gradFinalIter = 0;
	uint32_t Ndx = 1u, Ndy = 1u, Ndz = 0u, NiterAD = 1u, dimmu, inffi, Nlx = 1u, Nly = 1u, Nlz = 0u;
	bool med_no_norm = false, MBSREM_prepass = false, NLM_MRP = false, NLTV = false, NLRD = false, NLLange = false, NLM_anatomical = false, computeD = false,
		precondIm = false, precondMeas = false, computeM = false, RDPLargeNeighbor = false, UseL2Ball = true, NLLangeFiltered = false, filteringOrig = false, NLGGMRF = false, RDP_anatomical = false, 
		NLAdaptive = false;
	std::vector<bool> precondTypeMeas{ false, false }, precondTypeIm{ false, false, false, false, false, false, false };
	int64_t nProjections = 0LL;
	uint32_t nPriors = 0U, nMAP = 0U, nMAPML = 0U, nMLEM = 0U, nOS = 0U, nTot = 0U, nMAPOS = 0U, nPriorsTot = 0U, ng = 20U;
	std::vector<int32_t> mIt;
	std::vector<float> alphaCP, LCP, LCP2;
	float *rayShiftsDetector = nullptr, *rayShiftsSource = nullptr, *swivelAngles = nullptr;
	std::vector<float> alphaBB;
} Weighting;

// Struct for boolean operators indicating whether a certain method is selected
typedef struct RecMethods_ {
	bool OSEM = false, MRAMLA = false, RAMLA = false, ROSEM = false, RBI = false, DRAMA = false, COSEM = false, ECOSEM = false,
		ACOSEM = false, LSQR = false, CGLS = false, SART = false;
	bool MRP = false, Quad = false, Huber = false, L = false, FMH = false, WeightedMean = false, TV = false, AD = false, APLS = false, TGV = false, NLM = false, RDP = false, GGMRF = false, 
		ProxTV = false, ProxTGV = false, ProxRDP = false, ProxNLM = false, hyperbolic = false;
	bool OSLOSEM = false, MAPOSEM = false, MBSREM = false, BSREM = false, ROSEMMAP = false, ROSEMOSL = false, RBIMAP = false, RBIOSL = false,
		PKMA = false, FISTAL1 = false, SPS = false, PDHG = false, PDHGKL = false, FISTA = false, PDHGL1 = false, CV = false, PDDY = false, POCS = false, SAGA = false;
	bool MAP = false;
	bool CUSTOM = false;
	bool initAlg = false;
	bool CPType = false;
	bool OSL = false;
	bool FDK = false;
	bool BB = false;
	bool prior = false;
	uint32_t OSLCOSEM = 0u, MAPCOSEM = 0u;
} RecMethods;