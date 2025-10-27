// Structs to save buffer slots with Metal kernels and scalar inputs. This is a hybrid C++/MSL file with preprocessor directives. Possibly used with OpenCL/CUDA in the future for simpler code.

#if defined(__METAL_VERSION__) // MSL
#define FLOAT2_T float2
#define FLOAT3_T float3
#define UINT3_T uint3
#elif defined(__OPENCL_VERSION__) // OpenCL
#elif defined(__CUDACC__) // CUDA
#else // C++
#if defined(METAL) // C++ on MacOS
#pragma once
#include <cstdint>
#include <simd/simd.h>
#define FLOAT2_T simd::float2
#define FLOAT3_T simd::float3
#define UINT3_T simd::uint3
#elif defined(OPENCL) // C++ / OpenCL
#elif defined(CUDA) // C++ / CUDA
#endif
#endif

struct StaticScalarKernelParams { // Kernel scalar values that do not change with time step or (sub)iteration. See initializeKernel in ProjectorClass.h
    uint32_t nRowsD;
    uint32_t nColsD;
    float dPitchX;
    float dPitchY;
    float dL;
    float global_factor;
    float epps;
    uint32_t det_per_ring;
    float sigma_x;
	uint32_t rings;
    float coneOfResponseStdCoeffA;
    float coneOfResponseStdCoeffB;
    float coneOfResponseStdCoeffC;
    float tube_width;
    float cylRadiusProj3;
    float bmin;
    float bmax;
    float Vmax;
};

struct DynamicScalarKernelParams { // Kernel scalar values that do change with time step or (sub)iteration
    UINT3_T d_N;
    FLOAT3_T b;
    FLOAT2_T dSize;
    FLOAT3_T d;
    FLOAT3_T d_Scale;
    FLOAT3_T bmax;
    float orthWidth;
    long nProjections;
    unsigned char no_norm;
	unsigned long m_size;
	unsigned int currentSubset;
	int aa;
};

#undef FLOAT2_T
#undef FLOAT3_T
#undef UINT3_T