#pragma once
#ifdef CUDA
#include "cuda_error.h"
#define CTYPE3 float3
#define CTYPE2 float2
#elif defined(OPENCL)
#include "opencl_error.hpp"
#define CTYPE3 cl_float3
#define CTYPE2 cl_float2
#define ITYPE3 cl_int3
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
#include <cstdio>
#include <cstdint>
#include <fstream>
#ifdef MATLAB
#include "mexFunktio.h"
#endif

struct largeDimStruct {
	uint32_t NzOrig;
	float bzOrig;
	float bmaxZOrig;
	float d_Scale4ZOrig;
	int64_t imDimOrig;
	std::vector<uint32_t> Nz;
	std::vector<float> bz;
	std::vector<float> bmaxZ;
	std::vector<float> d_Scale4Z;
	std::vector<int64_t> imDim, cumDim;
};

typedef struct structForScalars {
	uint32_t projector_type = 1, attenuation_correction = 0, randoms_correction = 0, scatter = 0, normalization_correction = 0, 
		nColsD, nRowsD, size_z, subsets = 1, det_per_ring, Niter = 1, Nt = 1, subsetType = 0, nMultiVolumes = 0, nLayers = 1, 
		nRekos = 1, osa_iter0 = 0, nRekos2 = 0, subsetsUsed = 1, TOFsubsets = 1, Nxy = 0U, NxOrig = 0U, NyOrig = 0U, NzOrig = 0U, NxPrior = 0U, NyPrior = 0U, NzPrior = 0U,
		BPType = 1, FPType = 1, adaptiveType = 0, rings = 0, FISTAType = 0, maskFPZ = 1, maskBPZ = 1;
	uint32_t platform = 0;
	std::vector<uint32_t> Nx{ 1, 0, 0, 0, 0, 0, 0 }, Ny{ 1, 0, 0, 0, 0, 0, 0 }, Nz{ 1, 0, 0, 0, 0, 0, 0 };
	float crystal_size_z = 0.f, epps = 1e-6f, sigma_x = 0.f, tube_width = 0.f, bmin = 0.f, bmax = 0.f, Vmax = 0.f, global_factor = 1.f,
		dL = 0.f, flat = 0.f, cylRadiusProj3 = 0.f, DSC = 0.f;//, T = 0.f
	std::vector<float> dx{ 0.f, 0.f }, dy{ 0.f, 0.f }, dz{ 0.f, 0.f }, bx{ 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f }, by{ 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f }, 
		bz{ 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f };
	bool use_psf = false, TOF = false, SPECT = false, pitch = false, PET = false, meanFP = false, meanBP = false,
		maskFP = false, maskBP = false, orthXY = false, orthZ = false, CT = false, atomic_64bit = false, atomic_32bit = false, loadTOF = true, 
		saveIter = false, enforcePositivity = false, computeSensImag = false, useMAD = true, useImages = false, eFOV = false ,
		useExtendedFOV = false, use64BitIndices = false, TGV2D = false, multiResolution = false, offset = false, relaxScaling = false, 
		computeRelaxation = false, storeFP = false, deconvolution = false, CTAttenuation = true, largeDim = false, storeResidual = false, 
		useBuffers = true, useFDKWeights = false, indexBased = false, FISTAAcceleration = false;
	int64_t Nf = 0;
	std::vector<CTYPE3> d_Scale;
	std::vector<CTYPE3> d_Scale4;
	std::vector<CTYPE2> dSize;
	CTYPE2 dSizeBP;
	uint8_t raw = 0, fp = 0, listmode = 0;
	int8_t verbose = 0;
	uint16_t n_rays = 1, n_rays3D = 1;
	uint32_t g_dim_x = 0u, g_dim_y = 0u, g_dim_z = 0u, deblur_iterations = 0U;
	int64_t nBins = 1, nProjections = 0, numelY = 0, numelZ = 0, TOFSize = 0;
	std::vector<int64_t> im_dim{ 1, 0, 0, 0, 0, 0, 0 };
	size_t size_of_x, size_atten = 1, size_norm = 1, size_center_x, size_center_y, size_center_z, size_V = 1, size_scat = 1, koko = 0, sizeLOR, 
		sizeL, sizeXY, sizeZ, saveIterationsMiddle = 0ULL;
	uint32_t* saveNIter = nullptr;
	float* T = nullptr;
	float* V = nullptr;
	float* x_center = nullptr, *y_center = nullptr, *z_center = nullptr;
	float* gaussian = nullptr;
	float* TOFCenter = nullptr;
	uint64_t* pituus, length;
	std::vector<uint32_t> usedDevices;
	largeDimStruct lDimStruct;
	int64_t numMaskFP = 1, nProjectionsGlobal;
} scalarStruct;

#ifdef OPENCL
// Struct for the various estimates used in OpenCL kernels
typedef struct _OpenCL_im_vectors {
	cl::Buffer d_meanFP, d_meanBP;
	cl::Buffer	d_im;
	std::vector<cl::Buffer> d_rhs_os;
	cl::Image3D d_image_os, d_image_os_int;
} OpenCL_im_vectors;
#elif defined(CUDA)
typedef struct _CUDA_im_vectors {
	CUdeviceptr d_meanFP, d_meanBP;
	CUdeviceptr* d_im;
	std::vector<CUdeviceptr*> d_rhs_os;
	CUtexObject d_image_os, d_image_os_int;
} CUDA_im_vectors;
#else
struct CPUVectors {
	float* d_meanFP, *d_meanBP, *d_im_os;
	std::vector<float*> d_rhs_os;
};
#endif

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