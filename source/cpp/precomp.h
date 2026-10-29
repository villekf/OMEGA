#pragma once
#if defined(CUDA) || defined(HIP)
#include "cuda_error.h"
#define CTYPE3 float3
#define CTYPE2 float2
#elif defined(OPENCL)
#include "opencl_error.hpp"
#define CTYPE3 cl_float3
#define CTYPE2 cl_float2
#define ITYPE3 cl_int3
#elif defined(METAL)
#include "metal_error.hpp"
#include <simd/simd.h>
#define CTYPE3 simd::float3
#define CTYPE2 simd::float2
#define ITYPE3 simd::int3
#else
struct float3a {
	float x, y, z;
};
struct float2a {
	float x, y, z;
};
#define CTYPE3 float3a
#define CTYPE2 float2a
#endif
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdint>
#include <fstream>
#include <limits>
#include <cmath>
#ifdef MATLAB
#include "mexFunktio.h"
#endif
#ifdef METAL
#define NS_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#include "Metal.hpp"
#endif

struct largeDimStruct {
	uint32_t NzOrig;
	float bzOrig;
	float bmaxZOrig;
	float d_Scale4ZOrig;
	int64_t imDimOrig;
	std::vector<uint32_t> Nz, NzPr, startPr, endPr;
	std::vector<float> bz;
	std::vector<float> bmaxZ;
	std::vector<float> d_Scale4Z;
	std::vector<int64_t> imDim, cumDim, imDimPr, cumDimPr;
};

typedef struct structForScalars {
	uint32_t projector_type = 1, attenuation_correction = 0, randoms_correction = 0, scatter = 0, normalization_correction = 0, 
		nColsD, nRowsD, nHeads = 1, size_z, subsets = 1, det_per_ring, Niter = 1, Nt = 1, subsetType = 0, nMultiVolumes = 0, nLayers = 1,
		nRekos = 1, osa_iter0 = 0, timestep0 = 0, nRekos2 = 0, subsetsUsed = 1, timestepsUsed = 1, TOFsubsets = 1, Nxy = 0U, NxOrig = 0U, NyOrig = 0U, NzOrig = 0U, NxPrior = 0U, NyPrior = 0U, NzPrior = 0U,
		BPType = 1, FPType = 1, adaptiveType = 0, rings = 0, FISTAType = 0, maskFPZ = 1, normZ = 1, maskBPZ = 1, currentSubset = 0;
	uint32_t platform = 0;
	std::vector<uint32_t> Nx{ 1, 0, 0, 0, 0, 0, 0 }, Ny{ 1, 0, 0, 0, 0, 0, 0 }, Nz{ 1, 0, 0, 0, 0, 0, 0 };
	float crystal_size_z = 0.f, epps = 1e-6f, sigma_x = 0.f, tube_width = 0.f, bmin = 0.f, bmax = 0.f, Vmax = 0.f, global_factor = 1.f,
		dL = 0.f, flat = 0.f, cylRadiusProj3 = 0.f, DSC = 0.f, helicalRadius = 0.f;//, T = 0.f
	std::vector<float> dx{ 0.f, 0.f }, dy{ 0.f, 0.f }, dz{ 0.f, 0.f }, bx{ 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f }, by{ 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f }, 
		bz{ 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f };
    bool useHalf = false;
	bool use_psf = false, TOF = false, SPECT = false, pitch = false, PET = false, meanFP = false, meanBP = false,
		maskFP = false, maskBP = false, orthXY = false, orthZ = false, CT = false, atomic_64bit = false, atomic_32bit = false, loadTOF = true,
		saveIter = false, enforcePositivity = false, computeSensImag = false, useMAD = true, useImages = false, eFOV = false,
		useExtendedFOV = false, use64BitIndices = false, TGV2D = false, multiResolution = false, offset = false, relaxScaling = false,
		computeRelaxation = false, storeFP = false, deconvolution = false, CTAttenuation = true, largeDim = false, usingLinearizedData = false, storeResidual = false,
		useBuffers = true, useFDKWeights = false, indexBased = false, FISTAAcceleration = false, stochastic = false, useTotLength = true,
		useParallelBeam = false, useHelical = false;
	// Set this to false to explicitly disable fastPDHG
	// On MATLAB/Octave, you can do this just with options.fastPDHG = false;
	bool fastPDHG = true;
	// Build the projector type 2/3 kernels with the end point check of the tube (-DORTHCLIP), see orthClipNeeded below.
	// The default (true) is the conservative one, reconstructionAF()/reconstruction_multigpu() set it from the geometry
	bool orthClip = true;
	// True only when Nt > 1 AND some sensitivity-affecting input (image-based attenuation or the
	// multiplicative corrVector/extraCorr) genuinely differs per timestep. Computed once, early, in
	// reconstructionAF() and consumed there plus in functions.hpp (transferControl, backprojectionType6)
	// and subiterStep.h (computeOSEstimates) to pick vec.Summ's timestep slot: sIdx = sensPerTimestep ?
	// timestep : 0. When false (the common case, including all Nt == 1 runs), every timestep shares the
	// single sensitivity image computed from timestep 0 -- see reconstructionAF.h for why per-timestep
	// slots are otherwise unnecessary (normalization and measurement-domain attenuation are single
	// buffers shared by all timesteps).
	bool sensPerTimestep = false;
	int64_t Nf = 0;
	std::vector<CTYPE3> d_Scale;
	std::vector<CTYPE3> d_Scale4;
	std::vector<CTYPE2> dSize;
	CTYPE2 dSizeBP;
	uint8_t raw = 0, fp = 0, listmode = 0;
	int8_t verbose = 0;
	uint16_t n_rays = 1, n_rays3D = 1;
	uint32_t g_dim_x = 0u, g_dim_y = 0u, g_dim_z = 0u, deblur_iterations = 0U;
	// User-defined local/work-group (block) size per dimension. Any negative value means "use the
	// built-in default" for that dimension (see ProjectorClass::addProjector).
	int32_t localSize[3] = { -1, -1, -1 };
	// Compute the (spatial) regularization/prior only every regEveryIter-th time it is reached. A value
	// of 1 (or less) computes it every time (original behavior). The first and last time are always
	// computed. See subiterStep.h and iterStep.h.
	int32_t regEveryIter = 1;
	int64_t nBins = 1, nProjections = 0, numelY = 0, numelZ = 0, TOFSize = 0, seed = -1;
	std::vector<int64_t> im_dim{ 1 };
	size_t size_of_x, size_atten = 1, size_norm = 1, size_center_x, size_center_y, size_center_z, size_V = 1, size_scat = 1, kokoTOF = 0, kokoNonTOF = 0, sizeLOR,
		sizeL, sizeXY, sizeZ, saveIterationsMiddle = 0ULL;
	// Variables to check for correct size
	size_t size_maskFP = 0ULL, size_maskBP = 0ULL, size_meas = 0ULL;
	uint32_t* saveNIter = nullptr;
	float* T = nullptr;
	float* V = nullptr;
	float* x_center = nullptr, *y_center = nullptr, *z_center = nullptr;
	float* gaussian = nullptr;
	float* TOFCenter = nullptr;
	uint64_t* pituus, length;
	std::vector<uint32_t> usedDevices;
	largeDimStruct lDimStruct;
	float coneOfResponseStdCoeffA = 0.01f, coneOfResponseStdCoeffB = 0.01f, coneOfResponseStdCoeffC = 0.01f;
	float ellipseCenterX = 0.f, ellipseCenterY = 0.f, ellipseCenterZ = 0.f, ellipseRadiusX = 1.f, ellipseRadiusY = 1.f, ellipseRadiusZ = 1.f,
		ellipsePower = std::numeric_limits<float>::max(); // Box support, see extendRayToEllipse
} scalarStruct;

/// <summary>
/// Decides whether the projector types 2 and 3 are built with the end point check of the tube (-DORTHCLIP, see
/// orthogonalHelper3D in opencl_functions_orth3D.h). The tube of a ray ends at its end points, which only matters when a
/// voxel beyond an end point can be reached by the tube, i.e. when a detector/source is inside the image. Rule: the check is
/// enabled if and only if any detector/source coordinate point has its x and y inside the image box widened by m, and the z range of the
/// coordinates overlaps the box widened by m in z, where m is the tube reach (ORTH: tube_width, VOL: max(bmax, cylRadiusProj3))
/// plus the largest voxel dimension. The box is the union of the bounding boxes of all volumes (multi-resolution).
/// Reasoning: voxel centers are inside the box, so a voxel in the tube is closer than the reach to the line. For detectors
/// around the image (a ring) the line beyond an end point is outside the widened box (convex) and no voxel there is in the tube
/// With all coordinates outside the widened box the check is therefore skipped, which keeps the kernel (and its speed) unchanged.
/// The rule is conservative: it errs on the side of enabling the check. Non-finite coordinates (multi-layer padding) are ignored.
/// CT always enables it and SPECT never (the kernel does not use the check for SPECT).
/// </summary>
/// <param name="s various scalar parameters (geometry, volumes, projector types, tube parameters)"></param>
/// <param name="x detector coordinates (PET sinogram: 4 floats per bin, index-based/sensitivity: 2 floats per detector)"></param>
/// <param name="z_det z coordinates (PET sinogram: 2 floats per sinogram, 3 with multiple layers, otherwise one per ring)"></param>
/// <param name="listCoord list-mode event coordinates, 6 floats per event (non-index-based list-mode only)"></param>
/// <param name="pituus cumulative sum of the number of measurements per subset/timestep"></param>
inline bool orthClipNeeded(const scalarStruct& s, const float* x, const float* z_det, const float* listCoord, const int64_t* pituus) {
	if (s.SPECT)
		return false;
	if (s.CT)
		return true;
	// The reach of the tube
	float reach = 0.f;
	if (s.FPType == 2 || s.BPType == 2)
		reach = (std::max)(reach, s.tube_width);
	if (s.FPType == 3 || s.BPType == 3)
		reach = (std::max)(reach, (std::max)(s.bmax, s.cylRadiusProj3));
	// Union of the volumes
	const size_t nVol = (std::min)({ static_cast<size_t>(s.nMultiVolumes) + 1, s.Nx.size(), s.Ny.size(), s.Nz.size(), s.dx.size(), s.dy.size(), s.dz.size(),
		s.bx.size(), s.by.size(), s.bz.size() });
	if (nVol == 0ULL)
		return true;
	float boxLo[3] = { std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max() };
	float boxHi[3] = { std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest() };
	float voxel = 0.f;
	for (size_t ii = 0; ii < nVol; ii++) {
		const float b[3] = { s.bx[ii], s.by[ii], s.bz[ii] };
		const float d[3] = { s.dx[ii], s.dy[ii], s.dz[ii] };
		const float N[3] = { static_cast<float>(s.Nx[ii]), static_cast<float>(s.Ny[ii]), static_cast<float>(s.Nz[ii]) };
		for (int k = 0; k < 3; k++) {
			boxLo[k] = (std::min)(boxLo[k], b[k]);
			boxHi[k] = (std::max)(boxHi[k], b[k] + N[k] * d[k]);
			voxel = (std::max)(voxel, d[k]);
		}
	}
	const float m = reach + voxel;
	for (int k = 0; k < 3; k++) {
		boxLo[k] -= m;
		boxHi[k] += m;
	}
	// Coordinate points. xyInside: a point is inside the widened box in x and y, [zLo, zHi]: the z range of the coordinates
	bool anyXY = false, xyInside = false;
	float zLo = std::numeric_limits<float>::max(), zHi = std::numeric_limits<float>::lowest();
	auto addXY = [&](const float px, const float py) {
		if (!std::isfinite(px) || !std::isfinite(py))
			return;
		anyXY = true;
		if (px >= boxLo[0] && px <= boxHi[0] && py >= boxLo[1] && py <= boxHi[1])
			xyInside = true;
	};
	auto addZ = [&](const float pz) {
		if (!std::isfinite(pz))
			return;
		zLo = (std::min)(zLo, pz);
		zHi = (std::max)(zHi, pz);
	};
	auto decided = [&]() { return xyInside && zLo <= boxHi[2] && zHi >= boxLo[2]; };
	const bool listNonIndex = s.listmode > 0 && !s.indexBased;
	// Detector coordinates: xy pairs (the sinogram bins are two consecutive pairs), z per sinogram or per ring
	// In non-index-based list-mode these are the full detectors, which exist only for the sensitivity image
	if (x != nullptr && (!listNonIndex || s.computeSensImag))
		for (size_t k = 0; k + 1 < s.size_of_x; k += 2)
			addXY(x[k], x[k + 1]);
	if (z_det != nullptr && (!listNonIndex || s.computeSensImag)) {
		// PET sinograms: source and detector z of each sinogram, with multiple layers the third value is the layer combination
		const size_t stride = (s.listmode == 0 && !s.raw) ? (s.nLayers > 1 ? 3ULL : 2ULL) : 1ULL;
		const size_t cnt = stride == 1ULL ? 1ULL : 2ULL;
		for (size_t k = 0; k + cnt <= s.size_z; k += stride)
			for (size_t c = 0; c < cnt; c++)
				addZ(z_det[k + c]);
	}
	if (decided())
		return true;
	// List-mode events: the same range of events as is uploaded to the device
	if (listNonIndex && listCoord != nullptr && pituus != nullptr) {
		for (uint32_t tt = 0; tt < s.Nt; tt++) {
			for (uint32_t kk = s.osa_iter0; kk < s.subsetsUsed; kk++) {
				const size_t idx = static_cast<size_t>(kk) + static_cast<size_t>(tt) * static_cast<size_t>(s.subsets);
				for (int64_t ev = pituus[idx]; ev < pituus[idx + 1]; ev++) {
					const float* p = &listCoord[static_cast<size_t>(ev) * 6ULL];
					addXY(p[0], p[1]);
					addZ(p[2]);
					addXY(p[3], p[4]);
					addZ(p[5]);
					if (decided())
						return true;
				}
			}
		}
	}
	// No coordinates available: nothing can be excluded
	if (!anyXY)
		return true;
	// No z coordinates: the z range cannot exclude anything
	if (zLo > zHi)
		return xyInside;
	return decided();
}

#ifdef OPENCL
// Struct for the various estimates used in OpenCL kernels
typedef struct _OpenCL_im_vectors {
	cl::Buffer d_meanFP, d_meanBP;
	cl::Buffer	d_im;
	std::vector<cl::Buffer> d_rhs_os;
	cl::Image3D d_image_os, d_image_os_int;
} OpenCL_im_vectors;
#elif defined(CUDA) || defined(HIP)
typedef struct _CUDA_im_vectors {
	CUdeviceptr d_meanFP, d_meanBP;
	CUdeviceptr* d_im;
	std::vector<CUdeviceptr*> d_rhs_os;
	CUtexObject d_image_os, d_image_os_int;
} CUDA_im_vectors;
#elif defined(METAL)
typedef struct _METAL_im_vectors {
	NS::SharedPtr<MTL::Buffer> d_im;
    std::vector<NS::SharedPtr<MTL::Buffer>> d_rhs_os;
    NS::SharedPtr<MTL::Texture> d_image_os, d_image_os_int;
} METAL_im_vectors;
#else
struct CPUVectors {
	float* d_meanFP, *d_meanBP, *d_im_os;
	std::vector<float*> d_rhs_os;
};
#endif

// Optional device-resident inputs/outputs for implementation 5 (e.g. MATLAB gpuArray data)
// Every member is a device pointer valid in the backend's current context; nullptr means the
// corresponding host array is used instead, i.e. the original behavior
typedef struct _deviceIO {
	const void* im = nullptr;      // Forward projection input image (type 1)
	const void* meas = nullptr;    // Backprojection input measurements (type 2)
	void* output = nullptr;        // FP measurement output (type 1) or BP image output (type 2)
	void* sensIm = nullptr;        // BP sensitivity image output (type 2, only when no_norm == 0)
} deviceIO;

inline void mexPrint(const char* str) {
#ifdef MATLAB
	mexPrintf("%s\n", str);
	mexEvalString("pause(.0001);");
#else
	fprintf(stdout, "%s\n", str);
	fflush(stdout);
#endif
}

inline void mexPrintVar(const char* str, const int var) {
#ifdef MATLAB
	mexPrintf("%s%d\n", str, var);
	mexEvalString("pause(.0001);");
#else
	fprintf(stdout, "%s%d\n", str, var);
	fflush(stdout);
#endif
}

inline void mexPrintVarf(const char* str, const float var) {
#ifdef MATLAB
	mexPrintf("%s%f\n", str, var);
	mexEvalString("pause(.0001);");
#else
	fprintf(stdout, "%s%f\n", str, var);
	fflush(stdout);
#endif
}

template <typename T>
inline void mexPrintBase(const char* str, const T var) {
#ifdef MATLAB
	mexPrintf(str, var);
#else
	fprintf(stdout, str, var);
#endif
}

template <typename T, typename K>
inline void mexPrintBase(const char* str, const T var1, const K var2) {
#ifdef MATLAB
	mexPrintf(str, var1, var2);
#else
	fprintf(stdout, str, var1, var2);
#endif
}

template <typename T, typename K, typename C>
inline void mexPrintBase(const char* str, const T var1, const K var2, const C var3) {
#ifdef MATLAB
	mexPrintf(str, var1, var2, var3);
#else
	fprintf(stdout, str, var1, var2, var3);
#endif
}

inline void mexEval() {
#ifdef MATLAB
	mexEvalString("pause(.0001);");
#else
	fflush(stdout);
#endif
}

inline void mexWarning(const char* str) {
#ifdef MATLAB
	mexWarnMsgTxt(str);
#endif
}
