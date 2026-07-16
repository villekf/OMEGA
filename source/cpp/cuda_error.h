#pragma once
#if defined(HIP)
// ===================== HIP backend =====================
// HIP mirrors the CUDA driver / NVRTC API with renamed symbols. The reconstruction code is
// written against the CUDA driver API; this alias layer maps those names to their HIP/HIPRTC
// equivalents so the shared "#if defined(CUDA) || defined(HIP)" code paths compile under HIP.
#include <hip/hip_runtime.h>
#include <hip/hip_vector_types.h>
#include <hip/hiprtc.h>
#include <iostream>

// --- handle / scalar types ---
#define CUdeviceptr                                   hipDeviceptr_t
#define CUfunction                                    hipFunction_t
#define CUmodule                                      hipModule_t
#define CUresult                                      hipError_t
#define CUstream                                      hipStream_t
#define cuStreamCreate                                hipStreamCreateWithFlags
#define cuStreamDestroy                               hipStreamDestroy
#define cuStreamWaitEvent                             hipStreamWaitEvent
#define cuEventDestroy                                hipEventDestroy
#define CU_STREAM_NON_BLOCKING                        hipStreamNonBlocking
#define CU_EVENT_DISABLE_TIMING                       hipEventDisableTiming
// Runtime-API stream handle. Maps to the same type as the driver-API CUstream so that code using
// either name shares one stream type under HIP. NOTE: af/cuda.h must not also typedef cudaStream_t
// (see structs.h) or this alias would redefine hipStream_t.
#define cudaStream_t                                  hipStream_t
#define CUdevice                                      hipDevice_t
#define CUevent                                       hipEvent_t
#define CUtexObject                                   hipTextureObject_t
#define CUarray                                       hipArray_t
#define CUaddress_mode                                HIPaddress_mode
#define CUfilter_mode                                 HIPfilter_mode
#define CUarray_format                                hipArray_Format
#define CUmemorytype                                  hipMemoryType
#define CUresourcetype                                HIPresourcetype
#define CUresourceViewFormat                          HIPresourceViewFormat
#define CUDA_SUCCESS                                  hipSuccess
// --- descriptor structs ---
#define CUDA_ARRAY_DESCRIPTOR                         HIP_ARRAY_DESCRIPTOR
#define CUDA_ARRAY3D_DESCRIPTOR_st                    HIP_ARRAY3D_DESCRIPTOR
#define CUDA_MEMCPY2D                                 hip_Memcpy2D
#define CUDA_MEMCPY3D                                 HIP_MEMCPY3D
#define CUDA_RESOURCE_DESC                            HIP_RESOURCE_DESC
#define CUDA_RESOURCE_VIEW_DESC                       HIP_RESOURCE_VIEW_DESC
#define CUDA_TEXTURE_DESC                             HIP_TEXTURE_DESC
// --- driver functions ---
#define cuMemAlloc                                    hipMalloc
#define cuMemFree                                     hipFree
#define cuMemGetInfo                                  hipMemGetInfo
#define cuMemcpyHtoD                                  hipMemcpyHtoD
#define cuMemcpyDtoH                                  hipMemcpyDtoH
#define cuMemcpy2D                                    hipMemcpyParam2D
#define cuMemcpy3D                                    hipDrvMemcpy3D
#define cuMemcpy3DAsync                               hipDrvMemcpy3DAsync
#define cuModuleGetFunction                           hipModuleGetFunction
#define cuModuleLoadDataEx                            hipModuleLoadDataEx
#define cuModuleUnload                                hipModuleUnload
#define cuLaunchKernel                                hipModuleLaunchKernel
#define cuCtxSynchronize                              hipDeviceSynchronize
#define cuStreamSynchronize                           hipStreamSynchronize
#define cuDeviceGetAttribute                          hipDeviceGetAttribute
#define cuArrayCreate                                 hipArrayCreate
#define cuArray3DCreate                               hipArray3DCreate
#define cuArrayDestroy                                hipArrayDestroy
#define cuTexObjectCreate                             hipTexObjectCreate
#define cuTexObjectDestroy                            hipTexObjectDestroy
#define cuEventCreate                                 hipEventCreateWithFlags
#define cuEventRecord                                 hipEventRecord
#define cuEventSynchronize                            hipEventSynchronize
#define cuEventElapsedTime                            hipEventElapsedTime
// --- enums / flags ---
#define CU_EVENT_DEFAULT                              hipEventDefault
#define CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR  hipDeviceAttributeComputeCapabilityMajor
#define CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR  hipDeviceAttributeComputeCapabilityMinor
#define CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK hipDeviceAttributeMaxSharedMemoryPerBlock
#define CU_MEMORYTYPE_ARRAY                           hipMemoryTypeArray
#define CU_MEMORYTYPE_DEVICE                          hipMemoryTypeDevice
#define CU_MEMORYTYPE_HOST                            hipMemoryTypeHost
#define CU_RESOURCE_TYPE_ARRAY                        HIP_RESOURCE_TYPE_ARRAY
#define CU_RES_VIEW_FORMAT_FLOAT_1X32                 HIP_RES_VIEW_FORMAT_FLOAT_1X32
#define CU_RES_VIEW_FORMAT_UINT_1X8                   HIP_RES_VIEW_FORMAT_UINT_1X8
#define CU_TRSF_NORMALIZED_COORDINATES                HIP_TRSF_NORMALIZED_COORDINATES
#define CU_TRSF_READ_AS_INTEGER                       HIP_TRSF_READ_AS_INTEGER
#define CU_TR_ADDRESS_MODE_CLAMP                      HIP_TR_ADDRESS_MODE_CLAMP
#define CU_TR_FILTER_MODE_LINEAR                      HIP_TR_FILTER_MODE_LINEAR
#define CU_TR_FILTER_MODE_POINT                       HIP_TR_FILTER_MODE_POINT
#define CU_AD_FORMAT_FLOAT                            HIP_AD_FORMAT_FLOAT
#define CU_AD_FORMAT_UNSIGNED_INT8                    HIP_AD_FORMAT_UNSIGNED_INT8
// --- HIPRTC (runtime kernel compilation) ---
#define nvrtcProgram                                  hiprtcProgram
#define nvrtcResult                                   hiprtcResult
#define NVRTC_SUCCESS                                 HIPRTC_SUCCESS
#define NVRTC_ERROR_BUILTIN_OPERATION_FAILURE         HIPRTC_ERROR_BUILTIN_OPERATION_FAILURE
#define nvrtcCreateProgram                            hiprtcCreateProgram
#define nvrtcCompileProgram                           hiprtcCompileProgram
#define nvrtcDestroyProgram                           hiprtcDestroyProgram
#define nvrtcGetErrorString                           hiprtcGetErrorString
#define nvrtcGetProgramLog                            hiprtcGetProgramLog
#define nvrtcGetProgramLogSize                        hiprtcGetProgramLogSize
#define nvrtcGetPTX                                   hiprtcGetCode
#define nvrtcGetPTXSize                               hiprtcGetCodeSize

#define getErrorString(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(hipError_t code, const char* file, int line)
{
	if (code != hipSuccess)
		std::cerr << "GPUassert: " << hipGetErrorString(code) << ", " << file << ", line " << line << std::endl;
}
#else
// ===================== CUDA backend =====================
#include <nvrtc.h>
#include <vector_types.h>
#include <cuda.h>
#include <iostream>

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
#endif

#define TH 100000000000.f
#define TH32 100000.f
#define NVOXELS 8
#define NVOXELSHELICAL 1
#define NVOXELS5 1
#define NVOXELSFP 8
