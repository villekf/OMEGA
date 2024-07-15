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

//const char* getErrorString(CUresult error);

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