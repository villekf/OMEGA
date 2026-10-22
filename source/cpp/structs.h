#pragma once
// If true, displays a lot of debug messages
#define DEBUG false
// If true, uses the (correct) method to compute the inverse of measurement-based filtered PDHG, if false uses the same version as non-filtered one
// Only affects measurement-based filtering
#define FINVERSE true
// If false, uses the original EM preconditioner. If true, uses an EM preconditioner with the measurement data instead of ones.
#define SENSSCALE false
// If true, caches the NLM neighborhood in local memory
// Recommended to be true
#define FASTNLMLOCALCACHE true
// If true, outputs all the multi-resolution volumes (and the reconstruction itself) as different cell elements
#define CELL false
#include <algorithm>
#include <vector>
#include <cmath>
#include <numeric>
#include "precomp.h"
#ifdef AF
#include <arrayfire.h>
#if defined(CUDA)
#define AF_DEFINE_CUDA_TYPES 1
#include <af/cuda.h>
#include <cstring>
#elif defined(HIP)
// We must NOT define AF_DEFINE_CUDA_TYPES here: that makes af/cuda.h declare
// "typedef struct CUstream_st *cudaStream_t;", which our cudaStream_t -> hipStream_t alias
// (cuda_error.h) turns into a second, conflicting definition of hipStream_t (C2371). Instead we
// supply the only other type af/cuda.h needs (cublasMath_t) ourselves and let cudaStream_t
// resolve to hipStream_t through the alias, keeping a single consistent stream type.
typedef enum { CUBLAS_DEFAULT_MATH = 0, CUBLAS_TENSOR_OP_MATH = 1 } cublasMath_t;
#include <af/cuda.h>
#include <cstring>
#elif defined(OPENCL)
#include <af/opencl.h>
#elif defined(METAL)
#include <af/metal.h>
#endif
#endif

// Used in ProjectorClass.h (combined OpenCL and CUDA backends)
#if defined(OPENCL)
#define OCL_CHECK(STATUS, MSG, RETURN) do { \
    if ((STATUS) != CL_SUCCESS) { \
        getErrorString(STATUS); \
        mexPrint(MSG); \
        return RETURN; \
    } \
} while(0)
#elif defined(CUDA) || defined(HIP)
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
    // Arrays of size N_timesteps
    std::vector<af::array> dU; // gradient of spatial prior
    std::vector<af::array> dUt; // gradient of temporal prior
    std::vector<af::array> p0CP; 
    // Arrays of size N_timesteps x N_multivolumes
    std::vector<std::vector<af::array>> im_os; // estimate of image volume
    std::vector<std::vector<af::array>> rhs_os; // backprojection
    std::vector<std::vector<af::array>> uCP;
    std::vector<std::vector<af::array>> pCP;
    std::vector<std::vector<af::array>> wLSQR;
    std::vector<std::vector<af::array>> fLSQR;
    std::vector<std::vector<af::array>> rhsCP;
    std::vector<std::vector<af::array>> f0POCS;
    std::vector<std::vector<af::array>> uFISTA;
    std::vector<std::vector<af::array>> imBB;
    std::vector<std::vector<af::array>> gradBB;
    // Arrays of size N_timesteps x N_multivolumes x N_subsets
	std::vector<std::vector<std::vector<af::array>>> Summ; // sensitivity image
	// Never-seen voxel mask (see accumulateSeenMask in functions.hpp): [sensitivity slot][volume] boolean images that are true where at least one subset
	// has a sensitivity above the epps floor, and [slot][volume][subset] flags of the subsets that have been accumulated into it so far
	std::vector<std::vector<af::array>> seenMask;
	std::vector<std::vector<std::vector<char>>> seenDone;
    // Other arrays
    af::array C_co, rCGLS, meanFP, meanBP, apu, adapTypeA;
	std::vector<af::array> im_os_blurred, qProxTGV, vProxTGV, qProxTV, qProx, SAGASum;
	std::vector<af::array> fCGLS;
	std::vector<std::vector<af::array>> stochasticHelper;
} AF_im_vectors;
#endif

// Struct for various parameters, mainly various weights and coefficients
typedef struct Weighting_ {
#ifdef AF
	af::array tr_offsets, weights_quad, weights_TV, weights_huber, fmh_weights, a_L, weighted_weights, UU, Amin, weights_RDP;
	std::vector<af::array> preRef, gradF, RDPref;
    std::vector<std::vector<af::array>> dP, M, D;
	af::array gaussianNLM, gFilter, filter, filterIm, Ffilter;
	// Square root of the offset redundancy weights for each subset (nRowsD x 1 x number of projections), used for the forward projection
	// when inputScalars.offsetSqrt is enabled
	std::vector<af::array> offsetSqrtW;
	// The noise-weighted (NW-FBP) ramp filter bank, one filter per noise-weighting level (Nf x nLevelsNW)
	af::array filterNW;
	uint32_t nLevelsNW = 10U;
	af_flux_function FluxType;
	af_diffusion_eq DiffusionType;
#endif
	TVdata data;
	float *sigma_PKMA = nullptr, * uv = nullptr, * angles = nullptr, * NLM_ref = nullptr, * refImage = nullptr, * alphaPrecond = nullptr, * kerroin4 = nullptr, * weights = nullptr, *listCoord = nullptr, * RDP_ref = nullptr;
	uint32_t* rekot = nullptr;
    int32_t *distInt = nullptr, *distInt2 = nullptr;
	// Sub-voxel accurate versions of the above, used by the rotation-based SPECT projector.
	// distFloat is the distance from the collimator surface to the first image plane, in voxels,
	// while distFloat2 is the lateral shift of the detector panel, also in voxels
	float *distFloat = nullptr, *distFloat2 = nullptr;
	uint8_t* maskFP = nullptr, * maskBP = nullptr, * eFOVIndices = nullptr, *maskPrior = nullptr, *maskOffset = nullptr, *TOFIndices = nullptr;
	uint16_t* axIndex = nullptr, * trIndex = nullptr;
	uint32_t* detectorVector = nullptr;
	// Host arrays of the normalization and the multiplicative scatter/extra correction (corrVector), used by the SPECT projector type 6 (see type6Weights)
	const float* normHost = nullptr, * corrHost = nullptr;
	// The projector type 6 reads its per-view arrays (the geometry swivelAngles and blurPlanes*, the DetectorVector, the per-projection normalization) at the view of the
	// call plus the offset of the frame. frameStart6[tt] is the first view of the frame tt in an array that holds the views of all the frames one frame after the
	// other (frameStart6[Nt] = the total number of views), geomAllFrames6 and normAllFrames6 tell whether the geometry and the normalization are such arrays (when
	// they are for one frame, which every frame reuses, the offset is 0), see setFrameViewsType6. geomViews6 is the number of entries of the geometry arrays (0 = unknown)
	std::vector<int64_t> frameStart6;
	bool geomAllFrames6 = false, normAllFrames6 = false;
	int64_t geomViews6 = 0;
	// The backprojection masks of the SPECT projector type 6 of the volumes of a multi-resolution reconstruction, [frame][volume] (see prepareMaskBPType6)
	std::vector<std::vector<std::vector<uint8_t>>> maskBPVol;
	std::vector<float> epsilon_mramla, betaLSQR, alphaLSQR, thetaLSQR, rhoLSQR, phiLSQR, betaFISTA, tFISTA, tNFista;
    std::vector<std::vector<float>> tauCP, tauCP2, LCP, LCP2, alphaCP, sigmaCP, sigma2CP, thetaCP, lambda, lambdaFiltered, alphaM, alphaBB;
    float U = 1000000.f, h_ACOSEM = 1.f, TimeStepAD, KAD, w_sum = 0.f, h2 = 1.f, huber_delta = 0.f, ACOSEM_rhs = 0.f, h_ACOSEM_2 = 1.f, RDP_gamma = 1.f,
		dPitchX, dPitchY, gammaCGLS = 0.f, gammaTempCGLS = 0.f, alphaCGLS = 0.f, nuIEM = 0.f, alphaCPTV = 1.f,
		gradV1 = 0.f, gradV2 = 0.f, betaReg = 0.f, alpha0CPTGV = 1.f, alpha1CPTGV = 1.f, GGMRF_p = 0.f, GGMRF_q = 0.f, GGMRF_c = 0.f, GGMRF_pqc = 0.f,
		beta = 0.f, beta_temporal = 0.f, dtvg = 0.f, alphaPOCS = 0.2f, rMaxPOCS = 0.95f, POCSepps = 1e-4f, POCSalphaRed = 0.95f, NLAdaptiveConstant = 1e-5f, TemporalTVsmoothing = 1e-8f;
	uint32_t alku_fmh = 0u, mean_type = 0u, powerIterations = 0, derivType = 0, gradInitIter = 0, filterIter = 0, gradFinalIter = 0;
	uint32_t Ndx = 1u, Ndy = 1u, Ndz = 0u, NiterAD = 1u, dimmu, inffi, Nlx = 1u, Nly = 1u, Nlz = 0u;
	bool med_no_norm = false, MBSREM_prepass = false, NLM_MRP = false, NLTV = false, NLRD = false, NLLange = false, NLM_anatomical = false, computeD = false,
		precondIm = false, precondMeas = false, computeM = false, RDPLargeNeighbor = false, UseL2Ball = true, NLLangeFiltered = false, filteringOrig = false, NLGGMRF = false, RDP_anatomical = false, 
		NLAdaptive = false, NLMaxWeight = false, NLGM = false;
	std::vector<bool> precondTypeMeas{ false, false }, precondTypeIm{ false, false, false, false, false, false, false };
	int64_t nProjections = 0LL;
	uint32_t nPriors = 0U, nMAP = 0U, nMAPML = 0U, nMLEM = 0U, nOS = 0U, nTot = 0U, nMAPOS = 0U, nPriorsTot = 0U, ng = 20U;
	std::vector<int32_t> mIt;
	float *rayShiftsDetector = nullptr, *rayShiftsSource = nullptr, *swivelAngles = nullptr;
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
	// Noise-weighted FBP (Zeng 2013). When selected, FDK is set as well: NW-FBP follows the FDK
	// path exactly, except that the filtering is performed with the ray-by-ray noise-weighted filter bank
	bool NWFBP = false;
	bool BB = false;
	bool prior = false;
    bool TemporalSmoothness = false; // First derivative smoothness prior
    bool TemporalTV = false; // Total variations prior
	uint32_t OSLCOSEM = 0u, MAPCOSEM = 0u;
} RecMethods;
