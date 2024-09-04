#pragma once
#include <nvrtc.h>
#include <vector_types.h>
#include <cuda.h>
#include <iostream>

#define TH 100000000000.f
#define TH32 100000.f
#define NVOXELS 8
#define NVOXELS5 1
#define NVOXELSFP 8

// Source: https://stackoverflow.com/questions/14038589/what-is-the-canonical-way-to-check-for-errors-using-the-cuda-runtime-api
#define getErrorString(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(CUresult code, const char* file, int line)
{
	if (code != CUDA_SUCCESS)
	{
		const char* errstr;
		cuGetErrorString(code, &errstr);
		std::cerr << "GPUassert: " << errstr << ", " << file << ", line " << line << std::endl;
	}
}