#ifdef __METAL_VERSION__ // Compiling MSL
#define FLOAT2_kernelParams float2
#define FLOAT3_kernelParams float3
#define PTR_DEV device
#define PTR_CONST constant
#define UINT3_kernelParams uint3
#else // Compiling C++
#pragma once
#include <cstdint>
#include <simd/simd.h>
#define FLOAT2_kernelParams simd::float2
#define FLOAT3_kernelParams simd::float3
#define PTR_DEV
#define PTR_CONST
#define UINT3_kernelParams simd::uint3
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
    PTR_DEV float* d_rayShiftsDetector;
    PTR_DEV float* d_rayShiftsSource;
    float coneOfResponseStdCoeffA;
    float coneOfResponseStdCoeffB;
    float coneOfResponseStdCoeffC;
    float tube_width;
    float cylRadiusProj3;
    float bmin;
    float bmax;
    float Vmax;
    PTR_CONST float* d_TOFCenter;
    PTR_CONST float* d_V;
};

struct DynamicScalarKernelParams { // Kernel scalar values that do change with time step or (sub)iteration
    UINT3_kernelParams d_N;
    FLOAT3_kernelParams b;
    FLOAT2_kernelParams dSize;
    FLOAT3_kernelParams d;
    FLOAT3_kernelParams d_Scale;
    FLOAT3_kernelParams bmax;
    float orthWidth;
    long nProjections;
    unsigned char no_norm;
	unsigned long m_size;
	unsigned int currentSubset;
	int aa;
};