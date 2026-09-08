/******************************************************************************************************************************************
* Class object for forward and backward projections. Combined OpenCL, CUDA/HIP, and Metal version.
*
* Select the backend with the CUDA, HIP, METAL, or OPENCL preprocessor definition.
*
* Copyright (C) 2022-2026 Ville-Veikko Wettenhovi, Niilo Saarlemo
*
* This program is free software: you can redistribute it and/or modify  it under the terms of the GNU General Public License as published 
* by the Free Software Foundation, either version 3 of the License, or  (at your option) any later version.
*
* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License  along with this program. If not, see <https://www.gnu.org/licenses/>.
******************************************************************************************************************************************/
#pragma once
#if !defined(CUDA) && !defined(HIP) && !defined(METAL) && !defined(OPENCL)
#define OPENCL
#endif
#include "structs.h"
#include <array>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fstream>

// ======== OpenCL / CUDA / Metal compatibility aliases and macros ========
struct EmptyTextureArray {};
#if defined(CUDA) || defined(HIP)
using STATUS_t = CUresult;
using KERNELHANDLE_t = CUfunction;
using PROGRAMHANDLE_t = CUmodule;
using INT32_t = int;
using UINT32_t = unsigned int;
using WORKRANGE_t = std::array<UINT32_t, 3>;
using INT64_t = int64_t;
using UINT64_t = uint64_t;
using FLOAT2_t = float2;
using FLOAT3_t = float3;
using INT3_t = int3;
using UINT2_t = uint2;
using UCHAR_t = unsigned char;
using DEVBUFF_t = CUdeviceptr;
#if !defined(AF)
// Standalone implementation 5 uses CUDA driver pointers directly. ArrayFire
// builds retain the pointer-to-device-pointer representation it expects.
using AFDEVBUFF_t = CUdeviceptr;
#else
using AFDEVBUFF_t = CUdeviceptr*;
#endif
using TEX2D_t = CUtexObject;
using TEX3D_t = CUtexObject;
using TEXARRAY_t = CUarray;
#define SUCCESS_VALUE CUDA_SUCCESS
#elif defined(METAL)
using STATUS_t = int;
using KERNELHANDLE_t = NS::SharedPtr<MTL::ComputePipelineState>;
using PROGRAMHANDLE_t = NS::SharedPtr<MTL::Library>;
using INT32_t = int;
using UINT32_t = unsigned int;
using WORKRANGE_t = std::array<UINT32_t, 3>;
using INT64_t = int64_t;
using UINT64_t = uint64_t;
using FLOAT2_t = simd::float2;
using FLOAT3_t = simd::float3;
using INT3_t = simd::int3;
using UINT2_t = simd::uint2;
using UCHAR_t = unsigned char;
using DEVBUFF_t = NS::SharedPtr<MTL::Buffer>;
using AFDEVBUFF_t = NS::SharedPtr<MTL::Buffer>;
using TEX2D_t = NS::SharedPtr<MTL::Texture>;
using TEX3D_t = NS::SharedPtr<MTL::Texture>;
using TEXARRAY_t = EmptyTextureArray;
#define SUCCESS_VALUE 0
#elif defined(OPENCL)
using STATUS_t = cl_int;
using KERNELHANDLE_t = cl::Kernel;
using PROGRAMHANDLE_t = cl::Program;
using INT32_t = cl_int;
using UINT32_t = cl_uint;
using WORKRANGE_t = cl::NDRange;
using INT64_t = cl_long;
using UINT64_t = cl_ulong;
using FLOAT2_t = cl_float2;
using FLOAT3_t = cl_float3;
using INT3_t = cl_int3;
using UINT2_t = cl_uint2;
using UCHAR_t = cl_uchar;
using DEVBUFF_t = cl::Buffer;
using AFDEVBUFF_t = cl::Buffer;
using TEX2D_t = cl::Image2D;
using TEX3D_t = cl::Image3D;
using TEXARRAY_t = EmptyTextureArray;
#define SUCCESS_VALUE CL_SUCCESS
#endif
// Define one vector load format
#if defined(CUDA) || defined(HIP) || defined(METAL)
#define VEC_X(VEC) ((VEC).x)
#define VEC_Y(VEC) ((VEC).y)
#define VEC_Z(VEC) ((VEC).z)
#elif defined(OPENCL)
#define VEC_X(VEC) ((VEC).s[0])
#define VEC_Y(VEC) ((VEC).s[1])
#define VEC_Z(VEC) ((VEC).s[2])
#endif
// Kernel workgroup ranges, both 2D and 3D
#if defined(CUDA) || defined(HIP) || defined(METAL)
#define SET_RANGE2(RANGE, X, Y) do { \
	(RANGE)[0] = static_cast<UINT32_t>(X); \
	(RANGE)[1] = static_cast<UINT32_t>(Y); \
	(RANGE)[2] = 1U; \
} while(0)
#define SET_RANGE3(RANGE, X, Y, Z) do { \
	(RANGE)[0] = static_cast<UINT32_t>(X); \
	(RANGE)[1] = static_cast<UINT32_t>(Y); \
	(RANGE)[2] = static_cast<UINT32_t>(Z); \
} while(0)
#elif defined(OPENCL)
#define SET_RANGE2(RANGE, X, Y) do { \
	(RANGE) = cl::NDRange(static_cast<size_t>(X), static_cast<size_t>(Y)); \
} while(0)
#define SET_RANGE3(RANGE, X, Y, Z) do { \
	(RANGE) = cl::NDRange(static_cast<size_t>(X), static_cast<size_t>(Y), static_cast<size_t>(Z)); \
} while(0)
#endif
#if defined(CUDA) || defined(HIP)
#define SET_LAUNCH_RANGE3(RANGE, X, Y, Z, LOCAL_RANGE) \
	SET_RANGE3((RANGE), ((X) / (LOCAL_RANGE)[0]), ((Y) / (LOCAL_RANGE)[1]), ((Z) / (LOCAL_RANGE)[2]))
#elif defined(OPENCL) || defined(METAL)
#define SET_LAUNCH_RANGE3(RANGE, X, Y, Z, LOCAL_RANGE) SET_RANGE3((RANGE), (X), (Y), (Z))
#endif
#define SET_RANGE_Z(RANGE, Z) SET_RANGE3((RANGE), (RANGE)[0], (RANGE)[1], (Z))
// Get kernels
#if defined(CUDA) || defined(HIP)
#define GET_KERNEL(KVAR, PROG, NAME) status = cuModuleGetFunction(&KVAR, PROG, NAME)
#define KCHECK(MSG) CUDA_CHECK(status, MSG, status)
#elif defined(METAL)
#define GET_KERNEL(KVAR, PROG, NAME) do { \
	NS::Error* metalKernelError = nullptr; \
	NS::SharedPtr<MTL::Function> metalFunction; \
	if (!(PROG)) { \
		status = -1; \
		mexPrintBase("Metal library for function '%s' is empty\n", (NAME)); \
	} else { \
		metalFunction = NS::TransferPtr((PROG)->newFunction(NS::String::string((NAME), NS::ASCIIStringEncoding))); \
	} \
	if (!metalFunction) { \
		status = -1; \
		if (PROG) { \
			mexPrintBase("Metal function '%s' was not found in the compiled library\n", (NAME)); \
		} \
	} else { \
		(KVAR) = NS::TransferPtr(mtlDevice->newComputePipelineState(metalFunction.get(), &metalKernelError)); \
		status = (KVAR).get() ? SUCCESS_VALUE : -1; \
	} \
	if (status != SUCCESS_VALUE && metalFunction) { \
		const char* metalMsg = (metalKernelError && metalKernelError->localizedDescription()) \
			? metalKernelError->localizedDescription()->utf8String() : "unknown Metal pipeline error"; \
		mexPrintBase("Metal pipeline creation failed: %s\n", metalMsg); \
	} \
} while(0)
#define KCHECK(MSG) CHECK(status, MSG, status)
#elif defined(OPENCL)
#define GET_KERNEL(KVAR, PROG, NAME) KVAR = cl::Kernel(PROG, NAME, &status)
#define KCHECK(MSG) OCL_CHECK(status, MSG, -1)
#endif
// Get the kernel and check that it loaded
#define CREATE_KERNEL(KVAR, PROG, NAME, MSG) GET_KERNEL(KVAR, PROG, NAME); KCHECK(MSG)
// Unified status check (uses the backend success code); replaces paired CUDA_CHECK / OCL_CHECK
#define CHECK(STATUS, MSG, RETURN) do { if ((STATUS) != SUCCESS_VALUE) { getErrorString(STATUS); mexPrint(MSG); return RETURN; } } while(0)
// Allocate a device buffer. FLAGS is the OpenCL cl_mem_flags (ignored by CUDA, which has no
// context/flags). Host pointer is always NULL here (data is written separately)
// SIZE is the number of bytes to allocate in BUF
#if defined(CUDA) || defined(HIP)
#define ALLOC_BUFFER(BUF, FLAGS, SIZE) status = cuMemAlloc(&BUF, SIZE)
#elif defined(METAL)
#define ALLOC_BUFFER(BUF, FLAGS, SIZE) do { \
	(BUF) = NS::TransferPtr(mtlDevice->newBuffer(static_cast<NS::UInteger>(SIZE), (MTL::ResourceOptions)MTL::ResourceStorageModeShared)); \
	status = (BUF).get() ? SUCCESS_VALUE : -1; \
} while(0)
#elif defined(OPENCL)
#define ALLOC_BUFFER(BUF, FLAGS, SIZE) BUF = cl::Buffer(CLContext, FLAGS, SIZE, NULL, &status)
#endif
// Backend-neutral buffer access flags. CUDA and Metal ignore the OpenCL access mode
// because their allocation APIs do not encode it.
#if defined(CUDA) || defined(HIP) || defined(METAL)
#define BACKEND_BUFFER_READ_ONLY 0
#define BACKEND_BUFFER_READ_WRITE 0
#elif defined(OPENCL)
#define BACKEND_BUFFER_READ_ONLY CL_MEM_READ_ONLY
#define BACKEND_BUFFER_READ_WRITE CL_MEM_READ_WRITE
#endif
// Upload SIZE bytes from host SRC into device buffer BUF (non-blocking/asynchronous on OpenCL)
#if defined(CUDA) || defined(HIP)
#define WRITE_BUFFER(BUF, SIZE, SRC) status = cuMemcpyHtoD(BUF, SRC, SIZE)
#elif defined(METAL)
#define WRITE_BUFFER(BUF, SIZE, SRC) do { \
	if ((BUF).get() && (BUF)->contents()) { \
		std::memcpy((BUF)->contents(), (SRC), SIZE); \
		status = SUCCESS_VALUE; \
	} else { \
		status = -1; \
	} \
} while(0)
#elif defined(OPENCL)
#define WRITE_BUFFER(BUF, SIZE, SRC) status = CLCommandQueue[0].enqueueWriteBuffer(BUF, CL_FALSE, 0, SIZE, SRC)
#endif
// Download SIZE bytes from device buffer BUF into host destination DST
#if defined(CUDA) || defined(HIP)
#define READ_BUFFER(BUF, SIZE, DST) status = cuMemcpyDtoH((DST), (BUF), (SIZE))
#elif defined(METAL)
#define READ_BUFFER(BUF, SIZE, DST) do { \
	if ((BUF).get() && (BUF)->contents()) { \
		std::memcpy((DST), (BUF)->contents(), SIZE); \
		status = SUCCESS_VALUE; \
	} else { \
		status = -1; \
	} \
} while(0)
#elif defined(OPENCL)
#define READ_BUFFER(BUF, SIZE, DST) status = CLCommandQueue[0].enqueueReadBuffer(BUF, CL_FALSE, 0, SIZE, DST)
#endif
// Queue/stream synchronization
#if defined(CUDA) || defined(HIP)
#define FINISH_QUEUE(STATUS, MSG, RETURN) do { (STATUS) = cuCtxSynchronize(); CHECK((STATUS), MSG, RETURN); } while(0)
#elif defined(METAL)
#define FINISH_QUEUE(STATUS, MSG, RETURN) do { (STATUS) = SUCCESS_VALUE; } while(0)
#elif defined(OPENCL)
#define FINISH_QUEUE(STATUS, MSG, RETURN) do { (STATUS) = CLCommandQueue[0].finish(); CHECK((STATUS), MSG, RETURN); } while(0)
#endif
#if defined(CUDA) || defined(HIP)
using TimerPoint = CUevent;
#define INIT_TIMER(START, END) do { cuEventCreate(&(START), CU_EVENT_DEFAULT); cuEventCreate(&(END), CU_EVENT_DEFAULT); } while(0)
#define START_TIMER(START) cuEventRecord((START), CLCommandQueue[0])
#define STOP_TIMER(END) cuEventRecord((END), CLCommandQueue[0])
#define PRINT_TIMER(START, END, MSG) do { \
	cuEventSynchronize((END)); \
	float seconds = 0.f; \
	cuEventElapsedTime(&seconds, (START), (END)); \
	seconds /= 1000.f; \
	mexPrintBase(MSG, seconds); \
} while(0)
#elif defined(OPENCL) || defined(METAL)
using TimerPoint = std::chrono::steady_clock::time_point;
#define INIT_TIMER(START, END) do {} while(0)
#define START_TIMER(START) do { (START) = std::chrono::steady_clock::now(); } while(0)
#define STOP_TIMER(END) do { (END) = std::chrono::steady_clock::now(); } while(0)
#define PRINT_TIMER(START, END, MSG) do { \
	const std::chrono::duration<double> seconds = (END) - (START); \
	mexPrintBase(MSG, seconds.count()); \
} while(0)
#endif
#if defined(CUDA) || defined(HIP)
#define BACKEND_TEXTURE_POINT CUfilter_mode::CU_TR_FILTER_MODE_POINT
#define BACKEND_TEXTURE_LINEAR CUfilter_mode::CU_TR_FILTER_MODE_LINEAR
#define BACKEND_TEXTURE_DEFAULT_FLAGS 0
#define BACKEND_TEXTURE_READ_AS_INTEGER CU_TRSF_READ_AS_INTEGER
#define BACKEND_TEXTURE_NORMALIZED CU_TRSF_NORMALIZED_COORDINATES
#elif defined(METAL) ||  defined(OPENCL)
#define BACKEND_TEXTURE_POINT 0
#define BACKEND_TEXTURE_LINEAR 0
#define BACKEND_TEXTURE_DEFAULT_FLAGS 0
#define BACKEND_TEXTURE_READ_AS_INTEGER 0
#define BACKEND_TEXTURE_NORMALIZED 0
#endif
#if defined(CUDA) || defined(HIP)
#define RESIZE_TEXTURE_VECTOR(TEXTURES, ARRAYS, SIZE) do { (TEXTURES).resize(SIZE); (ARRAYS).resize(SIZE); } while(0)
#define RESIZE_TEXTURE_ARRAY(ARRAYS, SIZE) do { (ARRAYS).resize(SIZE); } while(0)
#define CREATE_FLOAT_TEXTURE3D_FROM_HOST(TEX, ARRAY, SRC, X_DIM, Y_DIM, DEPTH, FILTER, FLAGS) do { \
	const auto textureSpec = cudaFloatTextureSpec((Y_DIM), (X_DIM), (DEPTH), (FILTER), (FLAGS)); \
	status = createCudaTexture3DFromHost((TEX), (ARRAY), (SRC), textureSpec); \
} while(0)
#define CREATE_FLOAT_TEXTURE3D_FROM_DEVICE(TEX, ARRAY, SRC, X_DIM, Y_DIM, DEPTH, FILTER, FLAGS) do { \
	const auto textureSpec = cudaFloatTextureSpec((Y_DIM), (X_DIM), (DEPTH), (FILTER), (FLAGS)); \
	status = createCudaTexture3DFromDevice((TEX), (ARRAY), reinterpret_cast<CUdeviceptr>(SRC), textureSpec); \
} while(0)
#define CREATE_FLOAT_TEXTURE3D_EMPTY(TEX, ARRAY, WIDTH, HEIGHT, DEPTH) do { \
	status = SUCCESS_VALUE; \
} while(0)
// Mask call sites pass x/Nx first and y/Ny second; CUDA texture specs take height before width.
#define CREATE_MASK_TEXTURE2D_FROM_HOST(TEX, ARRAY, SRC, X_DIM, Y_DIM, FLAGS) do { \
	const auto textureSpec = cudaMaskTextureSpec((Y_DIM), (X_DIM), 1, (FLAGS)); \
	status = createCudaTexture2DFromHost((TEX), (ARRAY), (SRC), textureSpec); \
} while(0)
#define CREATE_MASK_TEXTURE3D_FROM_HOST(TEX, TEX3D, ARRAY, SRC, X_DIM, Y_DIM, VIEW_DEPTH, COPY_DEPTH, ARRAY_DEPTH, FLAGS) do { \
	auto textureSpec = cudaMaskTextureSpec((Y_DIM), (X_DIM), (VIEW_DEPTH), (FLAGS)); \
	textureSpec.copyDepth = (COPY_DEPTH); \
	textureSpec.arrayDepth = (ARRAY_DEPTH); \
	status = createCudaTexture3DFromHost((TEX), (ARRAY), (SRC), textureSpec); \
} while(0)
#elif defined(METAL)
#define RESIZE_TEXTURE_VECTOR(TEXTURES, ARRAYS, SIZE) do { (TEXTURES).resize(SIZE); } while(0)
#define RESIZE_TEXTURE_ARRAY(ARRAYS, SIZE) do { } while(0)
#define CREATE_FLOAT_TEXTURE3D_FROM_HOST(TEX, ARRAY, SRC, X_DIM, Y_DIM, DEPTH, FILTER, FLAGS) do { \
	(TEX) = createMetalFloatTextureFromHost((SRC), metalTextureSpec((X_DIM), (Y_DIM), (DEPTH), true)); \
	status = (TEX).get() ? SUCCESS_VALUE : -1; \
} while(0)
#define CREATE_FLOAT_TEXTURE3D_FROM_DEVICE(TEX, ARRAY, SRC, X_DIM, Y_DIM, DEPTH, FILTER, FLAGS) do { \
	const auto textureSpec = metalTextureSpec((X_DIM), (Y_DIM), (DEPTH), true); \
	if (!(SRC) || !(SRC)->contents()) { \
		status = -1; \
	} else { \
		if (!(TEX) || (TEX)->width() != textureSpec.width || (TEX)->height() != textureSpec.height || (TEX)->depth() != textureSpec.depth) \
			(TEX) = createMetalFloatTextureEmpty(textureSpec); \
		if (!(TEX)) { \
			status = -1; \
		} else { \
			const MTL::Region textureRegion(0, 0, 0, textureSpec.width, textureSpec.height, textureSpec.depth); \
			const NS::UInteger bytesPerRow = textureSpec.width * textureSpec.elementSize; \
			const NS::UInteger bytesPerImage = bytesPerRow * textureSpec.height; \
			(TEX)->replaceRegion(textureRegion, 0, 0, (SRC)->contents(), bytesPerRow, bytesPerImage); \
			status = SUCCESS_VALUE; \
		} \
	} \
} while(0)
#define CREATE_FLOAT_TEXTURE3D_EMPTY(TEX, ARRAY, WIDTH, HEIGHT, DEPTH) do { \
	(TEX) = createMetalFloatTextureEmpty(metalTextureSpec((WIDTH), (HEIGHT), (DEPTH), true)); \
	status = (TEX).get() ? SUCCESS_VALUE : -1; \
} while(0)
#define CREATE_MASK_TEXTURE2D_FROM_HOST(TEX, ARRAY, SRC, X_DIM, Y_DIM, FLAGS) do { \
	(TEX) = createMetalMaskTextureFromHost((SRC), metalTextureSpec((X_DIM), (Y_DIM), 1, false)); \
	status = (TEX).get() ? SUCCESS_VALUE : -1; \
} while(0)
#define CREATE_MASK_TEXTURE3D_FROM_HOST(TEX, TEX3D, ARRAY, SRC, X_DIM, Y_DIM, VIEW_DEPTH, COPY_DEPTH, ARRAY_DEPTH, FLAGS) do { \
	(TEX) = createMetalMaskTextureFromHost((SRC), metalTextureSpec((X_DIM), (Y_DIM), (VIEW_DEPTH), true)); \
	status = (TEX).get() ? SUCCESS_VALUE : -1; \
} while(0)
#elif defined(OPENCL)
#define RESIZE_TEXTURE_VECTOR(TEXTURES, ARRAYS, SIZE) do { (TEXTURES).resize(SIZE); } while(0)
#define RESIZE_TEXTURE_ARRAY(ARRAYS, SIZE) do { } while(0)
#define CREATE_FLOAT_TEXTURE3D_FROM_HOST(TEX, ARRAY, SRC, X_DIM, Y_DIM, DEPTH, FILTER, FLAGS) do { \
	(TEX) = TEX3D_t(CLContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, format, (X_DIM), (Y_DIM), (DEPTH), 0, 0, const_cast<void*>(static_cast<const void*>(SRC)), &status); \
} while(0)
#define CREATE_FLOAT_TEXTURE3D_FROM_DEVICE(TEX, ARRAY, SRC, X_DIM, Y_DIM, DEPTH, FILTER, FLAGS) do { \
	cl::detail::size_t_array textureRegion = { (X_DIM), (Y_DIM), (DEPTH) }; \
	(TEX) = TEX3D_t(CLContext, CL_MEM_READ_ONLY, format, (X_DIM), (Y_DIM), (DEPTH), 0, 0, NULL, &status); \
	OCL_CHECK(status, "Image creation failed\n", -1); \
	status = CLCommandQueue[0].enqueueCopyBufferToImage((SRC), (TEX), 0, origin, textureRegion); \
	OCL_CHECK(status, "Image copy failed\n", -1); \
} while(0)
#define CREATE_FLOAT_TEXTURE3D_EMPTY(TEX, ARRAY, WIDTH, HEIGHT, DEPTH) do { \
	(TEX) = TEX3D_t(CLContext, CL_MEM_READ_ONLY, format, (WIDTH), (HEIGHT), (DEPTH), 0, 0, NULL, &status); \
} while(0)
#define CREATE_MASK_TEXTURE2D_FROM_HOST(TEX, ARRAY, SRC, X_DIM, Y_DIM, FLAGS) do { \
	(TEX) = TEX2D_t(CLContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, formatMask, (X_DIM), (Y_DIM), 0, const_cast<void*>(static_cast<const void*>(SRC)), &status); \
} while(0)
#define CREATE_MASK_TEXTURE3D_FROM_HOST(TEX, TEX3D, ARRAY, SRC, X_DIM, Y_DIM, VIEW_DEPTH, COPY_DEPTH, ARRAY_DEPTH, FLAGS) do { \
	(TEX3D) = TEX3D_t(CLContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, formatMask, (X_DIM), (Y_DIM), (VIEW_DEPTH), 0, 0, const_cast<void*>(static_cast<const void*>(SRC)), &status); \
} while(0)
#endif
// Append one kernel argument VAR. CUDA pushes its address into the argument vector VEC; OpenCL
// sets it on KERNEL at the running index IDX, reporting any error inline via getErrorString.
#if defined(CUDA) || defined(HIP)
#define KARG(VEC, KERNEL, IDX, VAR) VEC.emplace_back(reinterpret_cast<void*>(&VAR))
#define KARG_SCALAR(VEC, KERNEL, IDX, VAR) KARG(VEC, KERNEL, IDX, VAR)
#define KARG_METAL_SLOT(IDX, SLOT) do {} while(0)
#elif defined(METAL)
#define KARG(VEC, KERNEL, IDX, VAR) metalSetKernelArg(encoder, (IDX), (VAR))
#define KARG_SCALAR(VEC, KERNEL, IDX, VAR) do {} while(0)
#define KARG_METAL_SLOT(IDX, SLOT) do { (IDX) = static_cast<UINT32_t>(SLOT); } while(0)
#elif defined(OPENCL)
#define KARG(VEC, KERNEL, IDX, VAR) getErrorString(KERNEL.setArg(IDX++, VAR))
#define KARG_SCALAR(VEC, KERNEL, IDX, VAR) KARG(VEC, KERNEL, IDX, VAR)
#define KARG_METAL_SLOT(IDX, SLOT) do {} while(0)
#endif
#if defined(CUDA) || defined(HIP) || defined(METAL)
#define ADD_OPT(OPTS, FLAG) OPTS.push_back(FLAG)
#elif defined(OPENCL)
#define ADD_OPT(OPTS, FLAG) OPTS += " " FLAG
#endif
// Build option selecting the forward projection. hipRTC force-includes hiprtc_runtime.h, which
// uses FP as a type name, so a command-line -DFP breaks that header. HIP passes -DOMEGA_FP
// instead; general_opencl_functions.h maps it back to FP after the hipRTC prelude.
#if defined(HIP)
#define FP_FLAG "-DOMEGA_FP"
#else
#define FP_FLAG "-DFP"
#endif
// Append an integer-valued build option "FLAG=VALUE". CUDA stores the formatted std::string
// directly in the (std::string) options vector; OpenCL appends it to the options string.
#if defined(CUDA) || defined(HIP) || defined(METAL)
#define ADD_OPT_INT(OPTS, FLAG, VALUE) OPTS.push_back(FLAG "=" + std::to_string(VALUE))
#elif defined(OPENCL)
#define ADD_OPT_INT(OPTS, FLAG, VALUE) OPTS += (" " FLAG "=" + std::to_string(VALUE))
#endif
#if defined(CUDA) || defined(HIP)
#define BACKEND_STR "CUDA"
#define BUILD_SUCCESS_VALUE NVRTC_SUCCESS
#elif defined(METAL)
#define BACKEND_STR "Metal"
#define BUILD_SUCCESS_VALUE SUCCESS_VALUE
#elif defined(OPENCL)
#define BACKEND_STR "OpenCL"
#define BUILD_SUCCESS_VALUE SUCCESS_VALUE
#endif
// ================================================================

#if defined(METAL)
#define SET_KERNEL_ARG_BYTES(kernel, bytes, size, index) (kernel)->setBytes((const void*)&(bytes), (NS::UInteger)(size), (index))

inline void metalSetKernelArg(const NS::SharedPtr<MTL::ComputeCommandEncoder>& encoder, UINT32_t& index, const DEVBUFF_t& buffer) {
	encoder->setBuffer(buffer.get(), 0, static_cast<NS::UInteger>(index++));
}

inline void metalSetKernelArg(const NS::SharedPtr<MTL::ComputeCommandEncoder>& encoder, UINT32_t& index, const TEX3D_t& texture) {
	encoder->setTexture(texture.get(), static_cast<NS::UInteger>(index++));
}

template <typename T>
inline void metalSetKernelArg(const NS::SharedPtr<MTL::ComputeCommandEncoder>& encoder, UINT32_t& index, const T& value) {
	encoder->setBytes(&value, static_cast<NS::UInteger>(sizeof(T)), static_cast<NS::UInteger>(index++));
}

inline NS::SharedPtr<NS::Dictionary> makeMetalPreprocessorMacros(const std::vector<std::string>& options) {
	std::vector<NS::Object*> keys;
	std::vector<NS::Object*> values;
	keys.reserve(options.size());
	values.reserve(options.size());
	for (const std::string& option : options) {
		if (option.rfind("-D", 0) != 0)
			continue;
		const std::string body = option.substr(2);
		const size_t equals = body.find('=');
		const std::string key = equals == std::string::npos ? body : body.substr(0, equals);
		if (key.empty())
			continue;
		keys.push_back(NS::String::string(key.c_str(), NS::UTF8StringEncoding));
		if (equals == std::string::npos)
			values.push_back(NS::Number::number(1));
		else {
			const std::string value = body.substr(equals + 1);
			char* endPtr = nullptr;
			const unsigned long long numericValue = std::strtoull(value.c_str(), &endPtr, 10);
			if (endPtr && *endPtr == '\0')
				values.push_back(NS::Number::number(numericValue));
			else
				values.push_back(NS::String::string(value.c_str(), NS::UTF8StringEncoding));
		}
	}
	return NS::RetainPtr(NS::Dictionary::dictionary(
		const_cast<NS::Object* const*>(values.data()),
		const_cast<NS::Object* const*>(keys.data()),
		static_cast<NS::UInteger>(keys.size())));
}

#endif

/// <summary>
/// Class object for forward and backward projections. Combined OpenCL, CUDA/HIP, and Metal version
/// </summary>
class ProjectorClass {
	//private:
#if defined(METAL)
	// Keep Metal C++ autoreleased objects alive until the projector releases its SharedPtr members.
	NS::SharedPtr<NS::AutoreleasePool> autoreleasePool;
#endif
		// Local size
	size_t local_size[3];
	size_t local_sizePrior[3];
#if defined(CUDA) || defined(HIP)
	// Contents of general_opencl_functions.h. For CUDA/HIP this is handed to NVRTC/HIPRTC as an actual
	// header, so the kernel sources #include it rather than having it concatenated in front of them (as is
	// still done for OpenCL).
	std::string nvrtcKernelHeader;
	// Contents of opencl_functions_orth3D.h, supplied as a distinct header. Non-empty (and passed to
	// NVRTC/HIPRTC) only for the forward/backward projector programs that use it; empty for the auxiliary
	// kernels so it is not included in that build.
	std::string nvrtcOrthHeader;
#endif // END CUDA
	// Kernel input indices
	UINT32_t kernelInd_MRAMLA = 0;
	UINT32_t kernelIndFP = 0;
	UINT32_t kernelIndBP = 0;
	UINT32_t kernelIndFPSubIter = 0;
	UINT32_t kernelIndBPSubIter = 0;
	UINT32_t kernelIndSens = 0;
	FLOAT3_t ellipseCenter, ellipseRadii;
	// Crystal pitch
	FLOAT2_t dPitch;
	// Image dimensions
	INT3_t d_NOrig, d_NPrior;
	// Values to add to the global size to make it divisible by local size
	size_t erotus[3];
	size_t erotusPrior[3];
	size_t erotusPriorEFOV[3];
	size_t erotusSens[3];
	// Local and global sizes
	WORKRANGE_t local, global, localPrior, globalPrior, globalPriorEFOV;
	struct ResourceState {
		bool useBuffers = true;
		bool xC = false;
		bool yC = false;
		bool zC = false;
		bool V = false;
		bool TOF = false;
		bool eFOV = false;
		bool GGMRF = false;
		bool maskFP = false;
		bool maskBP = false;
		bool atten = false;
		bool attenM = false;
		bool norm = false;
		bool extra = false;
		bool raw = false;
		bool subInd = false;
		bool priorMask = false;
		int NLMRef = 0;
		bool BPIm = false;
		bool proj5Im = false;
		bool auxMod = false;
		bool FPMod = false;
		bool BPMod = false;
		bool SensMod = false;
		bool xFull = false;
		bool zFull = false;
		bool offsetT = false;
		bool geom5 = false;
		bool indexBased = false;
		bool TOFIndex = false;
		bool angle = false;
		bool rayShifts = false;
		int zType = -1;
		int xSteps = -1;
		int zSteps = -1;
		int nSteps = 0;
		int eSteps = 0;
		int lSteps = 0;
		int iSteps = 0;
		int TOFSteps = 0;
		int aSteps = 0;
		int oSteps = 0;
		int g5Steps = 0;
		int tSteps = 1;
		int attenSize = 0;
	};
	ResourceState memAlloc;

	template <typename K, typename T>
	inline K make_vec3(T a, T b, T c) {
		K apu;
		VEC_X(apu) = a;
		VEC_Y(apu) = b;
		VEC_Z(apu) = c;
		return apu;
	}

#if defined(OPENCL)
	bool constantBuffer = false;

	// Get the OpenCL context for the current platform
	STATUS_t clGetPlatformsContext(const uint32_t platform, cl::Context & context, std::vector<cl::CommandQueue>&commandQueues, 
		const std::vector<uint32_t>&usedDevices, std::vector<cl::Device>&devices) {
		STATUS_t status = SUCCESS_VALUE;

		// Get the number of platforms 
		std::vector<cl::Platform> platforms;
		status = cl::Platform::get(&platforms);
		OCL_CHECK(status, "\n", status);
		if (DEBUG) {
			mexPrintBase("platforms.size() = %u\n", platforms.size());
			mexEval();
		}

		if (platforms.size() == 0) {
			std::cerr << "No platforms available!" << std::endl;
			status = -1;
			return status;
		}
		if (platform >= platforms.size()) {
			std::cerr << "The specified platform number is greater than the available platform numbers!" << std::endl;
			status = -1;
			return status;
		}
		if (DEBUG) {
			mexPrintBase("platform = %u\n", platform);
			mexEval();
		}

		// Get context properties from the chosen platform
		cl_context_properties properties[] = { CL_CONTEXT_PLATFORM, reinterpret_cast <cl_context_properties>(platforms[platform]()), 0 };

		// Create context from the chosen platform
		// If a single device was selected (options.cpu_to_gpu_factor = 0), use GPU if possible
		context = cl::Context(CL_DEVICE_TYPE_ALL, properties, NULL, NULL, &status);
		OCL_CHECK(status, "\n", status);
		// Get device IDs
		std::vector<cl::Device> devices2;
		status = context.getInfo(CL_CONTEXT_DEVICES, &devices2);
		OCL_CHECK(status, "\n", status);
		devices.push_back(devices2[usedDevices[0]]);
		if (DEBUG) {
			mexPrintBase("devices.size() = %u\n", devices.size());
			mexPrintBase("devices[0] = %u\n", devices[0]);
			mexEval();
		}

		// Create the command queues
		// Enable out of order execution (devices can compute kernels at the same time)
		for (size_t i = 0; i < devices.size(); i++) {
			//commandQueues.push_back(cl::CommandQueue(context, devices[i], CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, &status));
			commandQueues.push_back(cl::CommandQueue(context, devices[i], 0, &status));
			OCL_CHECK(status, "\n", status);
		}
		if (DEBUG) {
			mexPrintBase("commandQueues.size() = %u\n", commandQueues.size());
			mexEval();
		}

		for (UINT32_t i = 0; i < commandQueues.size(); i++) {
			commandQueues[i].finish();
		}

		return status;
	}
#endif // END CUDA

#if defined(CUDA) || defined(HIP)
	struct CudaTextureSpec {
		size_t height = 0;
		size_t width = 0;
		size_t arrayDepth = 1;
		size_t copyDepth = 1;
		size_t viewDepth = 1;
		size_t elementSize = sizeof(float);
		CUarray_format arrayFormat = CUarray_format::CU_AD_FORMAT_FLOAT;
		unsigned int channels = 1;
		CUresourceViewFormat viewFormat = CUresourceViewFormat::CU_RES_VIEW_FORMAT_FLOAT_1X32;
		CUfilter_mode filterMode = CUfilter_mode::CU_TR_FILTER_MODE_POINT;
		unsigned int flags = 0;
	};

	inline CudaTextureSpec cudaFloatTextureSpec(const size_t height, const size_t width, const size_t depth,
		const CUfilter_mode filterMode = CUfilter_mode::CU_TR_FILTER_MODE_POINT, const unsigned int flags = 0) const {
		CudaTextureSpec spec;
		spec.height = height;
		spec.width = width;
		spec.arrayDepth = depth;
		spec.copyDepth = depth;
		spec.viewDepth = depth;
		spec.filterMode = filterMode;
		spec.flags = flags;
		return spec;
	}

	inline CudaTextureSpec cudaMaskTextureSpec(const size_t height, const size_t width, const size_t depth,
		const unsigned int flags = CU_TRSF_READ_AS_INTEGER) const {
		CudaTextureSpec spec;
		spec.height = height;
		spec.width = width;
		spec.arrayDepth = depth;
		spec.copyDepth = depth;
		spec.viewDepth = depth;
		spec.elementSize = sizeof(uint8_t);
		spec.arrayFormat = CUarray_format::CU_AD_FORMAT_UNSIGNED_INT8;
		spec.viewFormat = CUresourceViewFormat::CU_RES_VIEW_FORMAT_UINT_1X8;
		spec.flags = flags;
		return spec;
	}

	inline STATUS_t createCudaTextureFromArray(TEX3D_t& texture, const TEXARRAY_t array, const CudaTextureSpec& spec, const bool is3D) const {
		CUDA_RESOURCE_DESC resDescLocal;
		CUDA_TEXTURE_DESC texDescLocal;
		CUDA_RESOURCE_VIEW_DESC viewDescLocal;
		std::memset(&resDescLocal, 0, sizeof(resDescLocal));
		std::memset(&texDescLocal, 0, sizeof(texDescLocal));
		std::memset(&viewDescLocal, 0, sizeof(viewDescLocal));

		resDescLocal.resType = CUresourcetype::CU_RESOURCE_TYPE_ARRAY;
		resDescLocal.res.array.hArray = array;
		texDescLocal.addressMode[0] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
		texDescLocal.addressMode[1] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
		if (is3D)
			texDescLocal.addressMode[2] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
		texDescLocal.filterMode = spec.filterMode;
		texDescLocal.flags = spec.flags;
		viewDescLocal.height = spec.height;
		viewDescLocal.width = spec.width;
		if (is3D)
			viewDescLocal.depth = spec.viewDepth;
		viewDescLocal.format = spec.viewFormat;
		return cuTexObjectCreate(&texture, &resDescLocal, &texDescLocal, &viewDescLocal);
	}

	inline STATUS_t createCudaTexture2DFromHost(TEX2D_t& texture, TEXARRAY_t& array, const void* source, const CudaTextureSpec& spec) const {
		CUDA_ARRAY_DESCRIPTOR arrDesc;
		std::memset(&arrDesc, 0, sizeof(arrDesc));
		arrDesc.Format = spec.arrayFormat;
		arrDesc.NumChannels = spec.channels;
		arrDesc.Height = spec.height;
		arrDesc.Width = spec.width;
		STATUS_t status = cuArrayCreate(&array, &arrDesc);
		if (status != CUDA_SUCCESS)
			return status;

		CUDA_MEMCPY2D copy;
		std::memset(&copy, 0, sizeof(copy));
		copy.srcMemoryType = CUmemorytype::CU_MEMORYTYPE_HOST;
		copy.srcHost = source;
		copy.srcPitch = spec.width * spec.elementSize;
		copy.dstMemoryType = CUmemorytype::CU_MEMORYTYPE_ARRAY;
		copy.dstArray = array;
		copy.WidthInBytes = spec.width * spec.elementSize;
		copy.Height = spec.height;
		status = cuMemcpy2D(&copy);
		if (status != CUDA_SUCCESS)
			return status;
		return createCudaTextureFromArray(texture, array, spec, false);
	}

	inline STATUS_t createCudaTexture3DFromHost(TEX3D_t& texture, TEXARRAY_t& array, const void* source, const CudaTextureSpec& spec) const {
		CUDA_ARRAY3D_DESCRIPTOR_st arrDesc;
		std::memset(&arrDesc, 0, sizeof(arrDesc));
		arrDesc.Format = spec.arrayFormat;
		arrDesc.NumChannels = spec.channels;
		arrDesc.Height = spec.height;
		arrDesc.Width = spec.width;
		arrDesc.Depth = spec.arrayDepth;
		STATUS_t status = cuArray3DCreate(&array, &arrDesc);
		if (status != CUDA_SUCCESS)
			return status;

		CUDA_MEMCPY3D copy;
		std::memset(&copy, 0, sizeof(copy));
		copy.srcMemoryType = CUmemorytype::CU_MEMORYTYPE_HOST;
		copy.srcHost = source;
		copy.srcPitch = spec.width * spec.elementSize;
		copy.srcHeight = spec.height;
		copy.dstMemoryType = CUmemorytype::CU_MEMORYTYPE_ARRAY;
		copy.dstArray = array;
		copy.WidthInBytes = spec.width * spec.elementSize;
		copy.Height = spec.height;
		copy.Depth = spec.copyDepth;
		status = cuMemcpy3D(&copy);
		if (status != CUDA_SUCCESS)
			return status;
		return createCudaTextureFromArray(texture, array, spec, true);
	}

	inline STATUS_t createCudaTexture3DFromDevice(TEX3D_t& texture, TEXARRAY_t& array, const CUdeviceptr source, const CudaTextureSpec& spec) const {
		CUDA_ARRAY3D_DESCRIPTOR_st arrDesc;
		std::memset(&arrDesc, 0, sizeof(arrDesc));
		arrDesc.Format = spec.arrayFormat;
		arrDesc.NumChannels = spec.channels;
		arrDesc.Height = spec.height;
		arrDesc.Width = spec.width;
		arrDesc.Depth = spec.arrayDepth;
		STATUS_t status = cuArray3DCreate(&array, &arrDesc);
		if (status != CUDA_SUCCESS)
			return status;

		CUDA_MEMCPY3D copy;
		std::memset(&copy, 0, sizeof(copy));
		copy.srcMemoryType = CUmemorytype::CU_MEMORYTYPE_DEVICE;
		copy.srcDevice = source;
		copy.srcPitch = spec.width * spec.elementSize;
		copy.srcHeight = spec.height;
		copy.dstMemoryType = CUmemorytype::CU_MEMORYTYPE_ARRAY;
		copy.dstArray = array;
		copy.WidthInBytes = spec.width * spec.elementSize;
		copy.Height = spec.height;
		copy.Depth = spec.copyDepth;
		status = cuMemcpy3D(&copy);
		if (status != CUDA_SUCCESS)
			return status;
		return createCudaTextureFromArray(texture, array, spec, true);
	}
#elif defined(METAL)
	struct MetalTextureSpec {
		NS::UInteger width = 0;
		NS::UInteger height = 0;
		NS::UInteger depth = 1;
		NS::UInteger elementSize = sizeof(float);
		MTL::PixelFormat pixelFormat = MTL::PixelFormat::PixelFormatR32Float;
		bool force3D = false;
	};

	inline MetalTextureSpec metalTextureSpec(const size_t width, const size_t height, const size_t depth, const bool force3D) const {
		MetalTextureSpec spec;
		spec.width = static_cast<NS::UInteger>(width);
		spec.height = static_cast<NS::UInteger>(height);
		spec.depth = static_cast<NS::UInteger>(depth);
		spec.force3D = force3D;
		return spec;
	}

	inline TEX3D_t createMetalTexture(const MetalTextureSpec& spec) const {
		if (!mtlDevice || spec.width == 0 || spec.height == 0 || spec.depth == 0)
			return nullptr;

		const bool is3D = spec.force3D || spec.depth > 1;
		// Reject unsupported dimensions before calling Metal. Without this check MATLAB would crash instead of returning cleanly.
		constexpr NS::UInteger max2DDimension = 16384;
		constexpr NS::UInteger max3DDimension = 2048;
		const bool invalidDimensions = is3D
			? (spec.width > max3DDimension || spec.height > max3DDimension || spec.depth > max3DDimension)
			: (spec.width > max2DDimension || spec.height > max2DDimension);
		if (invalidDimensions) {
			mexPrintBase("Requested Metal texture dimensions: %llu x %llu x %llu\n", static_cast<unsigned long long>(spec.width), static_cast<unsigned long long>(spec.height), static_cast<unsigned long long>(is3D ? spec.depth : 1));
			mexWarning("Metal 2D textures have a maximum dimension of 16384 and 3D textures have a maximum dimension of 2048. Reconstruction was stopped. Use buffers or reduce the size of the requested texture.");
			return nullptr;
		}
		NS::SharedPtr<MTL::TextureDescriptor> desc = NS::TransferPtr(MTL::TextureDescriptor::alloc()->init());
		desc->setTextureType(is3D ? MTL::TextureType::TextureType3D : MTL::TextureType::TextureType2D);
		desc->setPixelFormat(spec.pixelFormat);
		desc->setWidth(spec.width);
		desc->setHeight(spec.height);
		desc->setDepth(is3D ? spec.depth : 1);

		TEX3D_t texture =
			NS::TransferPtr(mtlDevice->newTexture(desc.get()));
		return texture;
	}

	inline TEX3D_t createMetalTextureFromHost(const void* source, const MetalTextureSpec& spec) const {
		if (!source)
			return nullptr;

		TEX3D_t texture = createMetalTexture(spec);
		if (!texture)
			return nullptr;

		const bool is3D = spec.force3D || spec.depth > 1;
		MTL::Region textureRegion(0, 0, 0, spec.width, spec.height, is3D ? spec.depth : 1);
		const NS::UInteger bytesPerRow = spec.width * spec.elementSize;
		const NS::UInteger bytesPerImage = bytesPerRow * spec.height;
		if (is3D)
			texture->replaceRegion(textureRegion, 0, 0, source, bytesPerRow, bytesPerImage);
		else
			texture->replaceRegion(textureRegion, 0, source, bytesPerRow);
		return texture;
	}

	inline TEX3D_t createMetalFloatTextureFromHost(const float* source, const MetalTextureSpec& spec) const {
		return createMetalTextureFromHost(source, spec);
	}

	inline TEX3D_t createMetalFloatTextureEmpty(const MetalTextureSpec& spec) const {
		return createMetalTexture(spec);
	}

	inline TEX3D_t createMetalFloatTextureFromBuffer(const DEVBUFF_t& source, const MetalTextureSpec& spec) const {
		if (!source || !source->contents())
			return nullptr;
		return createMetalTextureFromHost(source->contents(), spec);
	}

	inline TEX3D_t createMetalMaskTextureFromHost(const uint8_t* source, const MetalTextureSpec& spec) const {
		if (!source || spec.width == 0 || spec.height == 0 || spec.depth == 0)
			return nullptr;
		std::vector<float> mask(static_cast<size_t>(spec.width) * static_cast<size_t>(spec.height) * static_cast<size_t>(spec.depth));
		for (size_t ii = 0; ii < mask.size(); ++ii)
			mask[ii] = static_cast<float>(source[ii]);
		MetalTextureSpec floatSpec = spec;
		floatSpec.pixelFormat = MTL::PixelFormat::PixelFormatR32Float;
		floatSpec.elementSize = sizeof(float);
		return createMetalTextureFromHost(mask.data(), floatSpec);
	}

	inline int updateMetalImageTextureFromBuffer(const scalarStruct& inputScalars, const int ii) {
		if (!vec_opencl.d_im || !vec_opencl.d_im->contents()) {
			// Standalone projector calls upload image-mode input directly into
			// d_image_os. ArrayFire calls instead provide d_im and require the
			// cached texture to be refreshed below after every subset/volume.
			if (vec_opencl.d_image_os)
				return 0;
			mexPrint("Unable to create Metal image texture: missing input buffer and texture");
			return -1;
		}
		const MetalTextureSpec spec = metalTextureSpec(
			inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii], true);
		const size_t volume = static_cast<size_t>(ii);
		if (FPTexCache.size() <= volume) {
			FPTexCache.resize(volume + 1);
			imageCacheDims.resize((volume + 1) * 3, 0);
		}
		if (!FPTexCache[volume] ||
			imageCacheDims[volume * 3] != spec.width ||
			imageCacheDims[volume * 3 + 1] != spec.height ||
			imageCacheDims[volume * 3 + 2] != spec.depth) {
			FPTexCache[volume] = createMetalTexture(spec);
			imageCacheDims[volume * 3] = spec.width;
			imageCacheDims[volume * 3 + 1] = spec.height;
			imageCacheDims[volume * 3 + 2] = spec.depth;
		}
		if (!FPTexCache[volume]) {
			mexPrint("Unable to create Metal image texture");
			return -1;
		}
		const MTL::Region textureRegion(0, 0, 0, spec.width, spec.height, spec.depth);
		const NS::UInteger bytesPerRow = spec.width * spec.elementSize;
		const NS::UInteger bytesPerImage = bytesPerRow * spec.height;
		FPTexCache[volume]->replaceRegion(textureRegion, 0, 0,
			vec_opencl.d_im->contents(), bytesPerRow, bytesPerImage);
		vec_opencl.d_image_os = FPTexCache[volume];
		return 0;
	}
#endif // END CUDA/METAL texture helpers

	/// <summary>
	/// This function creates the backend programs for projection and auxiliary reconstruction kernels
	/// </summary>
#if defined(OPENCL)
	/// <param name="CLContext OpenCL context"></param>
	/// <param name="CLDeviceID OpenCL device ID"></param>
#endif // END CUDA
	/// <param name="programFP the program to store forward projection program"></param>
	/// <param name="programBP the program to store backprojection program"></param>
	/// <param name="programAux the program to store auxliary (such as priors) programs"></param>
	/// <param name="header_directory the location of the kernel and header files"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters"></param>
	/// <param name="MethodList specifies the algorithms and priors used"></param>
	/// <param name="w_vec specifies some of the special options used"></param>
	/// <param name="local_size the local size"></param>
	/// <returns></returns>
#if defined(CUDA) || defined(HIP)
	inline nvrtcResult createProgram(PROGRAMHANDLE_t & programFP, PROGRAMHANDLE_t & programBP,
		PROGRAMHANDLE_t & programAux,
#elif defined(METAL)
	inline STATUS_t createProgram(PROGRAMHANDLE_t & programFP, PROGRAMHANDLE_t & programBP,
		PROGRAMHANDLE_t & programAux, PROGRAMHANDLE_t & programSens,
#elif defined(OPENCL)
	inline STATUS_t createProgram(cl::Context & CLContext, cl::Device & CLDeviceID, cl::Program & programFP, cl::Program & programBP,
		cl::Program & programAux, cl::Program & programSens,
#endif // END CUDA
		const char* header_directory, scalarStruct & inputScalars, const RecMethods MethodList,
		const Weighting & w_vec, const size_t local_size[], const int type = -1) {

#if defined(CUDA) || defined(HIP)
		int compMajor = 0, compMinor = 0;
		cuDeviceGetAttribute(&compMajor, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, CUDeviceID[0]);
		cuDeviceGetAttribute(&compMinor, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, CUDeviceID[0]);
		nvrtcResult status = NVRTC_SUCCESS;
#elif defined(METAL)
		STATUS_t status = SUCCESS_VALUE;
#ifdef AF
		mtlDevice = NS::RetainPtr(afmtl::getDevice());
#else
		mtlDevice = NS::TransferPtr(MTL::CreateSystemDefaultDevice());
#endif
		if (!mtlDevice) {
			mexPrint("No Metal device available");
			return -1;
		}
#elif defined(OPENCL)
		STATUS_t status = SUCCESS_VALUE;
		std::string deviceName = CLDeviceID.getInfo<CL_DEVICE_VENDOR>(&status);
		std::string NV("NVIDIA Corporation");
		std::string AMD("Advanced Micro Devices, Inc.");
#endif // END CUDA
		std::string kernelFile = header_directory;
		std::string kernel_path, kernel_pathBP;
		std::string contentFP, contentBP;
		std::string contentAux;
#if defined(CUDA) || defined(HIP)
		std::vector<std::string> options;

#if defined(CUDA)
		std::string buffer0 = "--gpu-architecture=compute_" + std::to_string(compMajor) + std::to_string(compMinor);
		options.push_back(buffer0);
		options.push_back("-DCUDA");
#elif defined(HIP)
		// Compile the HIP kernels for the architecture of the current GPU (e.g. gfx90a).
		hipDeviceProp_t hipProp;
		hipGetDeviceProperties(&hipProp, CUDeviceID[0]);
		std::string buffer0 = std::string("--gpu-architecture=") + hipProp.gcnArchName;
		options.push_back(buffer0);
		options.push_back("-DHIP");
#endif
#elif defined(METAL)
		std::vector<std::string> options;
		ADD_OPT(options, "-DMETAL");
		if (inputScalars.atomic_32bit)
			ADD_OPT(options, "-DATOMIC32");
#elif defined(OPENCL)
		std::string options = "-cl-single-precision-constant";
		cl::string apu = CLDeviceID.getInfo<CL_DEVICE_EXTENSIONS>();
		cl::string apu2 = "cl_ext_float_atomics";
		ADD_OPT(options, "-DOPENCL");
#endif // END CUDA
		if (inputScalars.useMAD) {
#if defined(CUDA)
			options.push_back("--use_fast_math");
#elif defined(HIP)
			// HIPRTC is Clang-based and does not accept NVRTC's --use_fast_math; use the Clang flag instead.
			options.push_back("-ffast-math");
#elif defined(OPENCL)
			ADD_OPT(options, "-cl-fast-relaxed-math");
#endif // END CUDA
			ADD_OPT(options, "-DUSEMAD");
		}
		if ((inputScalars.useImages && inputScalars.FPType != 4 && inputScalars.FPType != 5 && inputScalars.BPType != 5) || 
			(inputScalars.FPType == 4 || inputScalars.FPType == 5 || inputScalars.BPType == 5)) {
			ADD_OPT(options, "-DUSEIMAGES");
			inputScalars.useBuffers = false;
			memAlloc.useBuffers = false;
		}
		if (inputScalars.useHalf)
			ADD_OPT(options, "-DHALF");
		std::ifstream sourceHeader(kernelFile + "general_opencl_functions.h");
		// Load the header text file
		std::string contentHeader((std::istreambuf_iterator<char>(sourceHeader)), std::istreambuf_iterator<char>());
#if defined(METAL)
		{
			const std::string kernelParamsFile = kernelFile.substr(0, kernelFile.size() - 7) + "cpp/kernelParams.hpp";
			std::ifstream sourceKernelParams(kernelParamsFile);
			std::string contentKernelParams((std::istreambuf_iterator<char>(sourceKernelParams)), std::istreambuf_iterator<char>());
			contentHeader = contentKernelParams + contentHeader;
		}
#endif // END METAL
		// Orthogonal/volume of intersection header, kept as a distinct header (see below). It is only used
		// by the forward/backward projector programs (projectorType123), never by the auxiliary kernels.
		std::string contentOrth;
		bool useOrth = false;
		// Load orthogonal/volume of intersection headers if applicable
		if (inputScalars.FPType == 2 || inputScalars.BPType == 2 || inputScalars.FPType == 3 || inputScalars.BPType == 3) {
			if (inputScalars.orthXY)
				ADD_OPT(options, "-DCRYSTXY");
			if (inputScalars.orthZ)
				ADD_OPT(options, "-DCRYSTZ");
			std::ifstream sourceHeader3(kernelFile + "opencl_functions_orth3D.h");
			std::string contentHeader3((std::istreambuf_iterator<char>(sourceHeader3)), std::istreambuf_iterator<char>());
			contentOrth = contentHeader3;
			useOrth = true;
#if defined(OPENCL)
		}
		if (NV.compare(deviceName) == 0)
			ADD_OPT(options, "-DNVIDIA");
		else if (AMD.compare(deviceName) == 0) {
			cl_bool is_integrated = CLDeviceID.getInfo<CL_DEVICE_HOST_UNIFIED_MEMORY>(&status);
			if (status == CL_SUCCESS && !is_integrated)
				ADD_OPT(options, "-DAMD");
		}
		else if (apu.find(apu2) != std::string::npos) {
			ADD_OPT(options, "-DINTEL");
#endif // END CUDA
		}

#if defined(CUDA) || defined(HIP)
		// Hand the headers to NVRTC/HIPRTC as actual header files; the kernel sources only #include them.
		nvrtcKernelHeader = contentHeader;
		nvrtcOrthHeader = useOrth ? contentOrth : std::string();
		std::string headerPrefix = "#include \"general_opencl_functions.h\"\n";
		if (useOrth)
			headerPrefix += "#include \"opencl_functions_orth3D.h\"\n";
		// The auxiliary kernels never use the orthogonal/volume header, so they only include the general one.
		const std::string headerPrefixAux = "#include \"general_opencl_functions.h\"\n";
#elif defined(OPENCL) || defined(METAL)
		const std::string headerPrefix = useOrth ? (contentHeader + contentOrth) : contentHeader;
		const std::string& headerPrefixAux = contentHeader;
#endif // END CUDA

		kernel_path = kernelFile;
		kernel_pathBP = kernelFile;
		if (inputScalars.FPType > 0 && inputScalars.FPType != 6) {
			if (inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
				kernel_path += "projectorType123.cl";
			}
			else if (inputScalars.FPType == 4)
				kernel_path += "projectorType4.cl";
			else if (inputScalars.FPType == 5)
				kernel_path += "projectorType5.cl";
			std::ifstream sourceFile(kernel_path.c_str());
			std::string contentFFP((std::istreambuf_iterator<char>(sourceFile)), std::istreambuf_iterator<char>());
			contentFP = headerPrefix + contentFFP;
		}
		if (inputScalars.BPType > 0 && inputScalars.BPType != 6) {
			if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
				kernel_pathBP += "projectorType123.cl";
			}
			else if (inputScalars.BPType == 4)
				kernel_pathBP += "projectorType4.cl";
			else if (inputScalars.BPType == 5)
				kernel_pathBP += "projectorType5.cl";
			std::ifstream sourceFileBP(kernel_pathBP.c_str());
			std::string contentFBP((std::istreambuf_iterator<char>(sourceFileBP)), std::istreambuf_iterator<char>());
			contentBP = headerPrefix + contentFBP;
		}
		if (inputScalars.useParallelBeam)
			ADD_OPT(options, "-DPARALLEL");
		else if (inputScalars.useHelical)
			ADD_OPT(options, "-DHELICAL");

		// Load the source text file
		// Set all preprocessor definitions
		const bool siddonVal = (inputScalars.FPType == 1 || inputScalars.BPType == 1 || inputScalars.FPType == 4 || inputScalars.BPType == 4) 
			? true : false;
#if defined(OPENCL)
		if (constantBuffer || (inputScalars.listmode > 0 && !inputScalars.indexBased))
			ADD_OPT(options, "-DUSEGLOBAL");
#elif defined(METAL)
		if (inputScalars.listmode > 0 && !inputScalars.indexBased)
			ADD_OPT(options, "-DUSEGLOBAL");
#endif // END CUDA
		if (inputScalars.raw == 1)
			ADD_OPT(options, "-DRAW");
		if (inputScalars.maskFP) {
			ADD_OPT(options, "-DMASKFP");
			if (inputScalars.maskFPZ > 1)
				ADD_OPT(options, "-DMASKFP3D");
		}
		if (inputScalars.SPECT && inputScalars.maskFPZ == inputScalars.nHeads &&
			(inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3 ||
			 inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3)) {
			ADD_OPT(options, "-DMASKFPBYDETECTOR");
		}
		if (inputScalars.normalization_correction && inputScalars.SPECT && inputScalars.normZ == inputScalars.nHeads &&
			(inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3 ||
			 inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3)) {
			ADD_OPT(options, "-DNORMBYDETECTOR");
		}
		if (inputScalars.maskBP) {
			ADD_OPT(options, "-DMASKBP");
			if (inputScalars.maskBPZ > 1)
				ADD_OPT(options, "-DMASKBP3D");
		}
		if (inputScalars.useTotLength)// && !inputScalars.SPECT)
			ADD_OPT(options, "-DTOTLENGTH");
		// For largeDim, projector type 4 needs the entire volume sizes/locations to correctly interpolate the ray
		if (inputScalars.largeDim)
			ADD_OPT(options, "-DLARGEDIM");
		if (inputScalars.CT && MethodList.FDK && inputScalars.useFDKWeights)
			ADD_OPT(options, "-DFDK");
		if (inputScalars.offset)
			ADD_OPT(options, "-DOFFSET");
		if (inputScalars.attenuation_correction == 1u && inputScalars.CTAttenuation)
			ADD_OPT(options, "-DATN");
		else if (inputScalars.attenuation_correction == 1u && !inputScalars.CTAttenuation)
			ADD_OPT(options, "-DATNM");
		if (inputScalars.normalization_correction == 1u)
			ADD_OPT(options, "-DNORM");
		if (inputScalars.scatter == 1u)
			ADD_OPT(options, "-DSCATTER");
		if (inputScalars.randoms_correction == 1u)
			ADD_OPT(options, "-DRANDOMS");
		if (inputScalars.nLayers > 1U) {
			if (inputScalars.listmode > 0 && inputScalars.indexBased)
				ADD_OPT_INT(options, "-DNLAYERS", inputScalars.nLayers);
			else
				ADD_OPT_INT(options, "-DNLAYERS", inputScalars.nProjections / (inputScalars.nLayers * inputScalars.nLayers));
		}
		if (inputScalars.TOF) {
			ADD_OPT(options, "-DTOF");
		}
		if (inputScalars.CT)
			ADD_OPT(options, "-DCT");
		else if (inputScalars.PET && inputScalars.listmode == 0)
			ADD_OPT(options, "-DPET");
		else if (inputScalars.SPECT) {
			ADD_OPT(options, "-DSPECT");
		}

		ADD_OPT_INT(options, "-DNBINS", inputScalars.nBins);
		if (inputScalars.listmode == 1)
			ADD_OPT(options, "-DLISTMODE");
		else if (inputScalars.listmode == 2)
			ADD_OPT(options, "-DLISTMODE2");
		if (inputScalars.listmode > 0 && inputScalars.indexBased)
			ADD_OPT(options, "-DINDEXBASED");
		if ((siddonVal && ((inputScalars.n_rays * inputScalars.n_rays3D) > 1)) || inputScalars.SPECT) {
			ADD_OPT_INT(options, "-DN_RAYS", inputScalars.n_rays * inputScalars.n_rays3D);
			ADD_OPT_INT(options, "-DN_RAYS2D", inputScalars.n_rays);
			ADD_OPT_INT(options, "-DN_RAYS3D", inputScalars.n_rays3D);
		}
		if (inputScalars.pitch)
			ADD_OPT(options, "-DPITCH");
		if (((inputScalars.subsets > 1 && (inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7))) 
			&& !inputScalars.CT && !inputScalars.SPECT && !inputScalars.PET && inputScalars.listmode == 0)
			ADD_OPT(options, "-DSUBSETS");
		if (local_size[1] > 0ULL) {
			ADD_OPT_INT(options, "-DLOCAL_SIZE", local_size[0]);
			ADD_OPT_INT(options, "-DLOCAL_SIZE2", local_size[1]);
		}
		else {
			ADD_OPT_INT(options, "-DLOCAL_SIZE", local_size[0]);
		}
		if (inputScalars.subsets > 1 && inputScalars.listmode == 0) {
			ADD_OPT_INT(options, "-DSTYPE", inputScalars.subsetType);
			ADD_OPT_INT(options, "-DNSUBSETS", inputScalars.subsets);
		}
		if (DEBUG) {
			mexPrintBase("path = %s\n", kernel_path.c_str());
			mexPrintBase("pathBP = %s\n", kernel_pathBP.c_str());
			mexPrintBase("file = %s\n", kernelFile.c_str());
			mexPrintBase("inputScalars.BPType = %u\n", inputScalars.BPType);
			mexPrintBase("inputScalars.FPType = %u\n", inputScalars.FPType);
			mexEval();
		}
		auto buildBackendProgram = [&](std::string& content, PROGRAMHANDLE_t& program, auto& buildOptions) {
#if defined(OPENCL)
			return buildProgram(inputScalars.verbose, content, CLContext, CLDeviceID, program, inputScalars.atomic_64bit, 
				inputScalars.atomic_32bit, buildOptions);
#else
			return buildProgram(inputScalars.verbose, content, program, buildOptions);
#endif // END CUDA
		};
		// Build projector program
		if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3 || inputScalars.FPType == 1 || 
			inputScalars.FPType == 2 || inputScalars.FPType == 3) {
#if defined(CUDA) || defined(HIP) || defined(METAL)
			std::vector<std::string> os_options = options;
#elif defined(OPENCL)
			std::string os_options = options;
#endif // END CUDA
			ADD_OPT(os_options, "-DSIDDON");
			ADD_OPT(os_options, "-DATOMICF");
#if defined(CUDA) || defined(HIP) || defined(METAL)
			std::vector<std::string> os_optionsFP = os_options;
#elif defined(OPENCL)
			std::string os_optionsFP = os_options;
#endif // END CUDA
			ADD_OPT(os_optionsFP, FP_FLAG);
			if (inputScalars.FPType == 3)
				ADD_OPT(os_optionsFP, "-DVOL");
			if (inputScalars.FPType == 2 || inputScalars.FPType == 3)
				ADD_OPT(os_optionsFP, "-DORTH");
			if (inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
				if (DEBUG) {
					mexPrint("Trying to build FP 1-3 program\n");
				}
				status = buildBackendProgram(contentFP, programFP, os_optionsFP);
				if (status == BUILD_SUCCESS_VALUE && DEBUG) {
					mexPrint("FP 1-3 program built\n");
				}
				else if (status != BUILD_SUCCESS_VALUE)
					return status;
				if (status == BUILD_SUCCESS_VALUE)
					memAlloc.FPMod = true;
				}
			if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
				if (DEBUG) {
					mexPrint("Trying to build BP 1-3 program\n");
				}
				ADD_OPT(os_options, "-DBP");
				if (inputScalars.BPType == 3)
					ADD_OPT(os_options, "-DVOL");
				if (inputScalars.BPType == 2 || inputScalars.BPType == 3)
					ADD_OPT(os_options, "-DORTH");
				status = buildBackendProgram(contentBP, programBP, os_options);
				if (status == BUILD_SUCCESS_VALUE && DEBUG) {
					mexPrint("BP 1-3 program built\n");
				}
				else if (status != BUILD_SUCCESS_VALUE)
					return status;
				if (status == BUILD_SUCCESS_VALUE)
					memAlloc.BPMod = true;
				}
			}
		if (inputScalars.FPType == 4 || inputScalars.BPType == 4) {
#if defined(CUDA) || defined(HIP) || defined(METAL)
			std::vector<std::string> os_options = options;
#elif defined(OPENCL)
			std::string os_options = options;
#endif // END CUDA
			if (inputScalars.FPType == 4)
				ADD_OPT(os_options, FP_FLAG);
			ADD_OPT(os_options, "-DPTYPE4");
			// Non-positive interpolation length selects the Joseph-type projector, where the ray is sampled once per voxel plane along
			// its dominant axis (the input dL is then unused)
			if (inputScalars.dL <= 0.f)
				ADD_OPT(os_options, "-DJOSEPH");
			if (!inputScalars.largeDim) {
				ADD_OPT_INT(os_options, "-DNVOXELS", NVOXELS);
				if (inputScalars.useHelical) {
					ADD_OPT_INT(os_options, "-DNVOXELSHELICAL", NVOXELSHELICAL);
				}
			}
			if (inputScalars.FPType == 4) {
				status = buildBackendProgram(contentFP, programFP, os_options);
				if (status == BUILD_SUCCESS_VALUE && DEBUG) {
					mexPrint("FP 4 program built\n");
				}
				else if (status != BUILD_SUCCESS_VALUE)
					return status;
				if (status == BUILD_SUCCESS_VALUE)
					memAlloc.FPMod = true;
			}
			if (inputScalars.BPType == 4) {
				os_options = options;
				if (inputScalars.CT) {
					ADD_OPT(os_options, "-DBP4");
					ADD_OPT_INT(os_options, "-DNVOXELS", NVOXELS);
#ifndef METAL
					// fastPDHG build options
					if (inputScalars.fastPDHG) {
						addFastPDHGOptions(os_options, inputScalars, w_vec, MethodList, inputScalars.useHelical ? NVOXELSHELICAL : NVOXELS);
					}
#endif
				}
				else {
					ADD_OPT(os_options, "-DPTYPE4");
					ADD_OPT(os_options, "-DATOMICF");
					if (inputScalars.dL <= 0.f)
						ADD_OPT(os_options, "-DJOSEPH");
				}
				ADD_OPT(os_options, "-DBP");
				status = buildBackendProgram(contentBP, programBP, os_options);
				if (status == BUILD_SUCCESS_VALUE && DEBUG) {
					mexPrint("BP 4 program built\n");
				}
				else if (status != BUILD_SUCCESS_VALUE)
					return status;
				if (status == BUILD_SUCCESS_VALUE)
					memAlloc.BPMod = true;
			}
		}
		if (inputScalars.FPType == 5 || inputScalars.BPType == 5) {
#if defined(CUDA) || defined(HIP) || defined(METAL)
			std::vector<std::string> os_options = options;
#elif defined(OPENCL)
			std::string os_options = options;
#endif // END CUDA
			ADD_OPT(os_options, "-DPROJ5");
			if (inputScalars.meanFP)
				ADD_OPT(os_options, "-DMEANDISTANCEFP");
			if (inputScalars.meanBP)
				ADD_OPT(os_options, "-DMEANDISTANCEBP");
			if (inputScalars.FPType == 5)
				ADD_OPT(os_options, FP_FLAG);
			if (inputScalars.BPType == 5)
				ADD_OPT(os_options, "-DBP");
			// BDD uses precomputed values, but you can also use it with the original on the fly calculations
			// If GEOM5 is not defined, the kernel uses the old format, otherwise it assumes precomputed values
#if !defined(METAL)
			if (inputScalars.BPType == 5 && inputScalars.CT && inputScalars.listmode == 0)
				ADD_OPT(os_options, "-DGEOM5");
			// Same as above, but for BPType 5
			if (inputScalars.BPType == 5 && inputScalars.CT && inputScalars.fastPDHG) {
				addFastPDHGOptions(os_options, inputScalars, w_vec, MethodList, inputScalars.pitch ? 1ULL : NVOXELS5);
			}
#endif
			if (inputScalars.pitch) {
				ADD_OPT_INT(os_options, "-DNVOXELS5", 1);
			}
			else {
				ADD_OPT_INT(os_options, "-DNVOXELS5", NVOXELS5);
			}
			ADD_OPT_INT(os_options, "-DNVOXELSFP", NVOXELSFP);
			if (inputScalars.FPType == 5) {
				status = buildBackendProgram(contentFP, programFP, os_options);
				if (status == BUILD_SUCCESS_VALUE && DEBUG) {
					mexPrint("FP 5 program built\n");
				}
				else if (status != BUILD_SUCCESS_VALUE)
					return status;
				if (status == BUILD_SUCCESS_VALUE)
					memAlloc.FPMod = true;
			}
			else {
				status = buildBackendProgram(contentBP, programBP, os_options);
				if (status == BUILD_SUCCESS_VALUE && DEBUG) {
					mexPrint("BP 5 program built\n");
				}
				else if (status != BUILD_SUCCESS_VALUE)
					return status;
				if (status == BUILD_SUCCESS_VALUE)
					memAlloc.BPMod = true;
			}
		}
		if (inputScalars.computeSensImag && inputScalars.listmode > 0) {
#if defined(CUDA) || defined(HIP) || defined(METAL)
			std::vector<std::string> os_options = options;
#elif defined(OPENCL)
			std::string os_options = options;
#endif // END CUDA
			ADD_OPT(os_options, "-DBP");
			ADD_OPT(os_options, "-DATOMICF");
			ADD_OPT(os_options, "-DSENS");
			if (inputScalars.BPType == 3)
				ADD_OPT(os_options, "-DVOL");
			if (inputScalars.BPType == 2 || inputScalars.BPType == 3)
				ADD_OPT(os_options, "-DORTH");
			if (inputScalars.BPType == 4) {
				ADD_OPT(os_options, "-DPTYPE4");
				ADD_OPT_INT(os_options, "-DNVOXELS", NVOXELS);
				if (inputScalars.dL <= 0.f)
					ADD_OPT(os_options, "-DJOSEPH");
			}
			else
				ADD_OPT(os_options, "-DSIDDON");
			status = buildBackendProgram(contentBP, programSens, os_options);
			if (status != BUILD_SUCCESS_VALUE)
				return status;
			if (status == BUILD_SUCCESS_VALUE)
				memAlloc.SensMod = true;
		}
		// Build prior programs
		if (MethodList.NLM || MethodList.MRP || MethodList.RDP || w_vec.precondTypeMeas[1] || w_vec.precondTypeIm[5]
			|| MethodList.TV || MethodList.APLS || MethodList.hyperbolic || MethodList.ProxTV || MethodList.ProxTGV || MethodList.PKMA || 
			MethodList.BSREM || MethodList.RAMLA || MethodList.MRAMLA || MethodList.MBSREM || MethodList.CPType || MethodList.ProxRDP || 
			MethodList.ProxNLM || MethodList.GGMRF || inputScalars.projector_type == 6 || type == 0) {
			if (DEBUG) {
				mexPrint("Building aux programs\n");
			}
#if defined(CUDA) || defined(HIP)
			std::vector<std::string> optionsAux;
			optionsAux.push_back(buffer0);
#if defined(HIP)
			optionsAux.push_back("-DHIP");
#else
			optionsAux.push_back("-DCUDA");
#endif
			if (inputScalars.useMAD) {
#if defined(HIP)
				// HIPRTC is Clang-based and does not accept NVRTC's --use_fast_math; use the Clang flag instead.
				optionsAux.push_back("-ffast-math");
#else
				optionsAux.push_back("--use_fast_math");
#endif
				optionsAux.push_back("-DUSEMAD");
			}
#elif defined(METAL)
			std::vector<std::string> optionsAux;
			ADD_OPT(optionsAux, "-DMETAL");
			if (inputScalars.useMAD)
				ADD_OPT(optionsAux, "-DUSEMAD");
#endif // END CUDA
			std::string auxKernelPath = kernelFile + "auxKernels.cl";
			std::ifstream sourceFileAux(auxKernelPath.c_str());
			std::string contentAAux((std::istreambuf_iterator<char>(sourceFileAux)), std::istreambuf_iterator<char>());
			// The auxiliary kernels only use the general header, not the orthogonal/volume header.
			contentAux = headerPrefixAux + contentAAux;
#if defined(CUDA) || defined(HIP)
			// Ensure the orthogonal/volume header is not supplied to (nor included in) the aux program.
			nvrtcOrthHeader.clear();
			if (inputScalars.use64BitIndices) {
				optionsAux.push_back("-DLTYPE=long long");
				optionsAux.push_back("-DLTYPE3=long3");
			}
#elif defined(METAL)
			if (inputScalars.use64BitIndices) {
				ADD_OPT(optionsAux, "-DLTYPE=long");
				ADD_OPT(optionsAux, "-DLTYPE3=long3");
			}
#elif defined(OPENCL)
			std::string optionsAux;
			optionsAux = "-cl-single-precision-constant";
			ADD_OPT(optionsAux, "-DOPENCL");
			if (inputScalars.useMAD) {
				ADD_OPT(optionsAux, "-cl-fast-relaxed-math");
				ADD_OPT(optionsAux, "-DUSEMAD");
			}
			if (inputScalars.use64BitIndices) {
				ADD_OPT(optionsAux, "-DLTYPE=long");
				ADD_OPT(optionsAux, "-DLTYPE3=long3");
			}
#endif // END CUDA
			if (inputScalars.largeDim)
				ADD_OPT(optionsAux, "-DLARGEDIM");
			if (inputScalars.useExtendedFOV)
				ADD_OPT(optionsAux, "-DEFOV");
			if (inputScalars.useImages)
				ADD_OPT(optionsAux, "-DUSEIMAGES");
			if (type == 2) {
				if (inputScalars.use_psf)
					ADD_OPT(optionsAux, "-DPSF");
			}
			else if (type == 0) {
				if (inputScalars.CT)
					ADD_OPT(optionsAux, "-DCT");
				if (inputScalars.randoms_correction)
					ADD_OPT(optionsAux, "-DRANDOMS");
				if (inputScalars.use_psf)
					ADD_OPT(optionsAux, "-DPSF");
			}
			else
				ADD_OPT(optionsAux, "-DAF");
			if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution)) {
				ADD_OPT(optionsAux, "-DMASKPRIOR");
				if (inputScalars.maskBPZ > 1)
					ADD_OPT(optionsAux, "-DMASKBP3D");
			}
			if (inputScalars.eFOV && !inputScalars.multiResolution)
				ADD_OPT(optionsAux, "-DEFOVZ");
			if (MethodList.MRP) {
				ADD_OPT(optionsAux, "-DMEDIAN");
				ADD_OPT_INT(optionsAux, "-DSEARCH_WINDOW_X", w_vec.Ndx);
				ADD_OPT_INT(optionsAux, "-DSEARCH_WINDOW_Y", w_vec.Ndy);
				ADD_OPT_INT(optionsAux, "-DSEARCH_WINDOW_Z", inputScalars.Nz[0] == 1 ? 0U : w_vec.Ndz);
			}
			if (MethodList.NLM) {
				ADD_OPT(optionsAux, "-DNLM_");
				if (w_vec.NLM_MRP) {
					ADD_OPT_INT(optionsAux, "-DNLTYPE", 2);
				}
				else if (w_vec.NLTV) {
					ADD_OPT_INT(optionsAux, "-DNLTYPE", 1);
				}
				else if (w_vec.NLRD) {
					ADD_OPT_INT(optionsAux, "-DNLTYPE", 3);
				}
				else if (w_vec.NLLange) {
					ADD_OPT_INT(optionsAux, "-DNLTYPE", 4);
				}
				else if (w_vec.NLLangeFiltered) {
					ADD_OPT_INT(optionsAux, "-DNLTYPE", 5);
				}
				else if (w_vec.NLGGMRF) {
					ADD_OPT_INT(optionsAux, "-DNLTYPE", 6);
				}
				else {
					ADD_OPT_INT(optionsAux, "-DNLTYPE", 0);
				}
				if (w_vec.NLAdaptive)
					ADD_OPT(optionsAux, "-DNLMADAPTIVE");
				if (w_vec.NLM_anatomical)
					ADD_OPT(optionsAux, "-DNLMREF");
				ADD_OPT_INT(optionsAux, "-DSWINDOWX", w_vec.Ndx);
				ADD_OPT_INT(optionsAux, "-DSWINDOWY", w_vec.Ndy);
				ADD_OPT_INT(optionsAux, "-DSWINDOWZ", w_vec.Ndz);
				ADD_OPT_INT(optionsAux, "-DPWINDOWX", w_vec.Nlx);
				ADD_OPT_INT(optionsAux, "-DPWINDOWY", w_vec.Nly);
				ADD_OPT_INT(optionsAux, "-DPWINDOWZ", w_vec.Nlz);
			}
			if (MethodList.GGMRF) {
				ADD_OPT(optionsAux, "-DGGMRF");
				ADD_OPT_INT(optionsAux, "-DSWPRIORTYPE", 1);
				ADD_OPT_INT(optionsAux, "-DSWINDOWX", w_vec.Ndx);
				ADD_OPT_INT(optionsAux, "-DSWINDOWY", w_vec.Ndy);
				ADD_OPT_INT(optionsAux, "-DSWINDOWZ", w_vec.Ndz);
			}
			if (MethodList.hyperbolic) {
				ADD_OPT(optionsAux, "-DHYPER");
				ADD_OPT_INT(optionsAux, "-DSWPRIORTYPE", 2);
				ADD_OPT_INT(optionsAux, "-DSWINDOWX", w_vec.Ndx);
				ADD_OPT_INT(optionsAux, "-DSWINDOWY", w_vec.Ndy);
				ADD_OPT_INT(optionsAux, "-DSWINDOWZ", w_vec.Ndz);
			}
			if (MethodList.RDP) {
				ADD_OPT(optionsAux, "-DRDP");
				if (w_vec.RDPLargeNeighbor) {
					ADD_OPT(optionsAux, "-DRDPCORNERS");
					ADD_OPT_INT(optionsAux, "-DSWINDOWX", w_vec.Ndx);
					ADD_OPT_INT(optionsAux, "-DSWINDOWY", w_vec.Ndy);
					ADD_OPT_INT(optionsAux, "-DSWINDOWZ", w_vec.Ndz);
				}
				if (w_vec.RDP_anatomical)
					ADD_OPT(optionsAux, "-DRDPREF");
			}
			if (MethodList.ProxRDP && w_vec.RDPLargeNeighbor)
				ADD_OPT(optionsAux, "-DRDPCORNERS");
			if (MethodList.TV && !w_vec.data.TV_use_anatomical) {
				ADD_OPT(optionsAux, "-DTVGRAD");
				if (w_vec.data.TVtype == 6)
					ADD_OPT(optionsAux, "-DTVW1");
				else if (w_vec.data.TVtype == 4)
					ADD_OPT(optionsAux, "-DSATV");
				else if (w_vec.data.TVtype == 2)
					ADD_OPT(optionsAux, "-DJPTV");
				if (w_vec.derivType > 0) {
					ADD_OPT_INT(optionsAux, "-DDIFFTYPE", w_vec.derivType);
				}
			}
			else if ((MethodList.TV && w_vec.data.TV_use_anatomical) || MethodList.APLS) {
				ADD_OPT(optionsAux, "-DTVGRAD");
				if (w_vec.data.TVtype == 1)
					ADD_OPT(optionsAux, "-DANATOMICAL1");
				else if (w_vec.data.TVtype == 2)
					ADD_OPT(optionsAux, "-DANATOMICAL2");
				else if (w_vec.data.TVtype == 5 || MethodList.APLS)
					ADD_OPT(optionsAux, "-DANATOMICAL3");
				if (w_vec.derivType > 0) {
					ADD_OPT_INT(optionsAux, "-DDIFFTYPE", w_vec.derivType);
				}
			}
			if (MethodList.ProxTV) {
				ADD_OPT(optionsAux, "-DPROXTV");
				if (w_vec.UseL2Ball)
					ADD_OPT(optionsAux, "-DL2");
				if (w_vec.derivType > 0) {
					ADD_OPT_INT(optionsAux, "-DDIFFTYPE", w_vec.derivType);
				}
			}
			if (MethodList.ProxTGV || MethodList.TGV) {
				ADD_OPT(optionsAux, "-DPROXTV");
				ADD_OPT(optionsAux, "-DPROXTGV");
				if (w_vec.UseL2Ball)
					ADD_OPT(optionsAux, "-DL2");
				if (w_vec.derivType > 0) {
					ADD_OPT_INT(optionsAux, "-DDIFFTYPE", w_vec.derivType);
				}
				if (!inputScalars.TGV2D)
					ADD_OPT(optionsAux, "-DTGVZ");
			}
			if (MethodList.ProxRDP)
				ADD_OPT(optionsAux, "-DPROXRDP");
			if (local_sizePrior[1] > 0ULL) {
				ADD_OPT_INT(optionsAux, "-DLOCAL_SIZE", local_sizePrior[0]);
				ADD_OPT_INT(optionsAux, "-DLOCAL_SIZE2", local_sizePrior[1]);
				ADD_OPT_INT(optionsAux, "-DLOCAL_SIZE3", local_sizePrior[2]);
			}
			else {
				ADD_OPT_INT(optionsAux, "-DLOCAL_SIZE", local_sizePrior[0]);
			}
			if (MethodList.PKMA)
				ADD_OPT(optionsAux, "-DPKMA");
			else if (MethodList.MBSREM || MethodList.MRAMLA)
				ADD_OPT(optionsAux, "-DMBSREM");
			else if (MethodList.BSREM || MethodList.RAMLA)
				ADD_OPT(optionsAux, "-DBSREM");
			else if (MethodList.CPType) {
				ADD_OPT(optionsAux, "-DPDHG");
				if (inputScalars.subsets > 1)
					ADD_OPT(optionsAux, "-DSUBSETS");
			}
			if (inputScalars.projector_type == 6)
				ADD_OPT(optionsAux, "-DROTATE");
#if defined(CUDA) || defined(HIP)
			status = buildProgram(inputScalars.verbose, contentAux, programAux, optionsAux);
			if (status == NVRTC_SUCCESS && DEBUG) {
				mexPrint("Aux program built\n");
			}
			else if (status != NVRTC_SUCCESS)
				return status;
			if (status == NVRTC_SUCCESS)
				memAlloc.auxMod = true;
#elif defined(METAL)
			if (MethodList.CPType || inputScalars.projector_type == 6) {
				status = buildProgram(inputScalars.verbose, contentAux, programAux, optionsAux);
				if (status != SUCCESS_VALUE)
					return status;
				memAlloc.auxMod = true;
			}
			else {
				programAux = nullptr;
				status = SUCCESS_VALUE;
			}
#elif defined(OPENCL)
			status = buildProgram(inputScalars.verbose, contentAux, CLContext, CLDeviceID, programAux, inputScalars.atomic_64bit, 
				inputScalars.atomic_32bit, optionsAux);
#endif // END CUDA
		}
		if (DEBUG) {
			mexPrintBase("status = %u\n", status);
			mexPrintBase("w_vec.NLM_MRP = %u\n", w_vec.NLM_MRP);
			mexPrintBase("w_vec.NLTV = %u\n", w_vec.NLTV);
			mexPrintBase("w_vec.NLRD = %u\n", w_vec.NLRD);
			mexPrintBase("w_vec.NLLange = %u\n", w_vec.NLLange);
			mexPrintBase("w_vec.NLLangeFiltered = %u\n", w_vec.NLLangeFiltered);
			mexEval();
		}
		return status;
			}

	/// <summary>
	/// Builds one backend program from the supplied source
	/// </summary>
	/// <param name="verbose the level of verbosity"></param>
	/// <param name="contentFP program code"></param>
#if defined(OPENCL)
	/// <param name="CLContext OpenCL context"></param>
	/// <param name="CLDeviceID OpenCL device ID"></param>
#endif // END CUDA
	/// <param name="program the program where to store the built program"></param>
#if defined(OPENCL)
	/// <param name="atomic_64bit are 64-bit (int64) atomics used"></param>
	/// <param name="atomic_32bit are 32-bit (int) atomics used"></param>
#endif // END CUDA
	/// <param name="options preprocessor values for the build"></param>
	/// <returns></returns>
#if defined(METAL)
	inline STATUS_t buildProgram(const int8_t verbose, std::string& content, PROGRAMHANDLE_t& program, std::vector<std::string>& options) {
		if (!mtlDevice) {
			mexPrint("No Metal device available");
			return -1;
		}
		if (DEBUG || verbose >= 3) {
			for (size_t ll = 0; ll < options.size(); ll++)
				mexPrintBase("%s ", options[ll].c_str());
			mexPrintBase("%s\n", "");
		}
		NS::Error* err = nullptr;
		NS::SharedPtr<MTL::CompileOptions> compileOptions = NS::TransferPtr(MTL::CompileOptions::alloc()->init());
		compileOptions->setMathMode(std::find(options.begin(), options.end(), "-DUSEMAD") != options.end() ? MTL::MathModeFast : MTL::MathModeRelaxed);
		NS::SharedPtr<NS::Dictionary> macros = makeMetalPreprocessorMacros(options);
		compileOptions->setPreprocessorMacros(macros.get());
		program = NS::TransferPtr(mtlDevice->newLibrary(
			NS::String::string(content.c_str(), NS::UTF8StringEncoding),
			compileOptions.get(),
			&err));
		if (!program) {
			const char* msg = (err && err->localizedDescription())
				? err->localizedDescription()->utf8String()
				: "unknown Metal compile error";
			mexPrintBase("newLibrary failed: %s\n", msg);
			return -1;
		}
		if (DEBUG || verbose >= 3)
			mexPrint("Metal program built\n");
		return SUCCESS_VALUE;
	}
#else
#if defined(CUDA) || defined(HIP)
	inline nvrtcResult buildProgram(const int8_t verbose, std::string & content, CUmodule & module, std::vector<std::string>&options) {
		nvrtcResult status = NVRTC_SUCCESS;
		STATUS_t status2 = SUCCESS_VALUE;
		nvrtcProgram program;
		if (DEBUG || verbose >= 3) {
			for (int ll = 0; ll < options.size(); ll++)
				mexPrintBase("%s ", options[ll].c_str());
			mexPrintBase("%s\n", "");
#elif defined(OPENCL)
	inline STATUS_t buildProgram(const int8_t verbose, std::string contentFP, cl::Context & CLContext, cl::Device & CLDeviceID, 
		cl::Program & program, bool& atomic_64bit, const bool atomic_32bit, std::string options) {
		STATUS_t status = SUCCESS_VALUE;
		size_t pituus;
		if (atomic_64bit) {
			pituus = options.length();
			options += " -DCAST=long";
			options += " -DATOMIC";
			ADD_OPT_INT(options, "-DTH", TH);
#endif // END CUDA
		}
#if defined(CUDA) || defined(HIP)
		// Supply the headers as actual headers (resolved by the #include directives in the source) instead
		// of concatenating them into the source string. opencl_functions_orth3D.h is only passed when it is
		// used (i.e. for the projector programs; it is empty for the auxiliary kernels).
		std::vector<const char*> headerSources, headerNames;
		headerSources.push_back(nvrtcKernelHeader.c_str());
		headerNames.push_back("general_opencl_functions.h");
		if (!nvrtcOrthHeader.empty()) {
			headerSources.push_back(nvrtcOrthHeader.c_str());
			headerNames.push_back("opencl_functions_orth3D.h");
		}
		status = nvrtcCreateProgram(&program, content.c_str(), "32bit", static_cast<int>(headerNames.size()), headerSources.data(), headerNames.data());
		if (status != NVRTC_SUCCESS) {
			std::cerr << nvrtcGetErrorString(status) << std::endl;
			return status;
#elif defined(OPENCL)
		else if (atomic_32bit) {
			options += " -DCAST=int";
			options += " -DATOMIC32";
			ADD_OPT_INT(options, "-DTH", TH32);
#endif // END CUDA
		}
#if defined(CUDA) || defined(HIP)
		// Build the program
		std::vector<const char*> optionsC;
		optionsC.reserve(options.size());
		for (const std::string& opt : options)
			optionsC.push_back(opt.c_str());
		status = nvrtcCompileProgram(program, optionsC.size(), optionsC.data());
		// Build log in case of failure
		if (status != NVRTC_SUCCESS) {
			std::cerr << nvrtcGetErrorString(status) << std::endl;
			mexPrint("Failed to build CUDA program. Build log: \n");
			size_t len;
			char* buffer;
			nvrtcGetProgramLogSize(program, &len);
			buffer = (char*)calloc(len, sizeof(size_t));
			nvrtcGetProgramLog(program, buffer);
			mexPrintBase("%s\n", buffer);
			free(buffer);
			nvrtcDestroyProgram(&program);
			return status;
#elif defined(OPENCL)
		else {
			options += " -DCAST=float";
#endif // END CUDA
		}
#if defined(CUDA) || defined(HIP)
		else if (verbose > 1)
			mexPrint("CUDA program built\n");
			size_t ptxSize;
			status = nvrtcGetPTXSize(program, &ptxSize);
			if (status != NVRTC_SUCCESS) {
				std::cerr << nvrtcGetErrorString(status) << std::endl;
				return status;
#elif defined(OPENCL)
		if (DEBUG || verbose >= 3)
			mexPrintBase("%s\n", options.c_str());
		if (atomic_64bit) {
			cl::string apu = CLDeviceID.getInfo<CL_DEVICE_EXTENSIONS>();
			cl::string apu2 = "cl_khr_int64_base_atomics";
			// NOTE: find returns npos when the extension is missing, size_t is unsigned (var < 0 is never true)
			size_t var = apu.find(apu2);
			if (var == cl::string::npos) {
				options.erase(pituus, options.size() + 1);
				options += " -DCAST=float";
				status = -1;
#endif // END CUDA
			}
#if defined(CUDA) || defined(HIP)
			char* ptx = new char[ptxSize];
			status = nvrtcGetPTX(program, ptx);
			if (status != NVRTC_SUCCESS) {
				std::cerr << nvrtcGetErrorString(status) << std::endl;
				return status;
#elif defined(OPENCL)
			else {
				std::vector<std::string> testi;
				testi.push_back(contentFP);
				cl::Program::Sources source(testi);
				program = cl::Program(CLContext, source);
				status = program.build({ CLDeviceID }, options.c_str());
				if (status == CL_SUCCESS && (DEBUG || verbose >= 3)) {
					mexPrint("OpenCL program (64-bit atomics) built\n");
#endif // END CUDA
				}
#if defined(CUDA) || defined(HIP)
				status2 = cuModuleLoadDataEx(&module, ptx, 0, 0, 0);
				CUDA_CHECK(status2, "\n", NVRTC_ERROR_BUILTIN_OPERATION_FAILURE);
#elif defined(OPENCL)
				else if (status != CL_SUCCESS) {
					mexPrint("Failed to build 64-bit atomics program.\n");
#endif // END CUDA
					if (DEBUG) {
#if defined(CUDA) || defined(HIP)
						mexPrintBase("ptxSize = %u\n", ptxSize);
#elif defined(OPENCL)
						getErrorString(status);
						std::vector<cl::Device> dev;
						CLContext.getInfo(CL_CONTEXT_DEVICES, &dev);
						for (int ll = 0; ll < dev.size(); ll++) {
							cl_build_status status = program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(dev[ll]);
							if (status != CL_BUILD_ERROR)
								continue;
							std::string name = dev[ll].getInfo<CL_DEVICE_NAME>();
							std::string buildlog = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(dev[ll]);
							mexPrintBase("Build log for %s:\n %s", name.c_str(), buildlog.c_str());
#endif // END CUDA
						}
#if defined(CUDA) || defined(HIP)
						// Destroy the program.
						status = nvrtcDestroyProgram(&program);
						if (status != NVRTC_SUCCESS) {
							std::cerr << nvrtcGetErrorString(status) << std::endl;
							return status;
#endif // END CUDA
						}
#if defined(CUDA) || defined(HIP)
						delete[] ptx;
#elif defined(OPENCL)
						options.erase(pituus, options.size() + 1);
						options += " -DCAST=float";
					}
				}
			}
		else
			status = -1;
		// If not, use 32-bit atomic add (float)
		if (status != CL_SUCCESS) {
			status = CL_SUCCESS;
			atomic_64bit = false;
			std::vector<std::string> testi;
			testi.push_back(contentFP);
			cl::Program::Sources source(testi);
			program = cl::Program(CLContext, source);
			status = program.build({ CLDeviceID }, options.c_str());
			if (status == CL_SUCCESS && (DEBUG || verbose >= 3)) {
				mexPrint("OpenCL program built\n");
			}
			else if (status != CL_SUCCESS) {
				mexPrint("Failed to build OpenCL program.\n");
				getErrorString(status);
				std::vector<cl::Device> dev;
				CLContext.getInfo(CL_CONTEXT_DEVICES, &dev);
				for (int ll = 0; ll < dev.size(); ll++) {
					cl_build_status status = program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(dev[ll]);
					if (status != CL_BUILD_ERROR)
						continue;
					std::string name = dev[ll].getInfo<CL_DEVICE_NAME>();
					std::string buildlog = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(dev[ll]);
					mexPrintBase("Build log for %s:\n %s", name.c_str(), buildlog.c_str());
				}
			}
		}
#endif // END CUDA
		return status;
	}
#endif // END METAL buildProgram

		/// <summary>
	/// Creates the necessary backend kernels from the input programs
	/// </summary>
	/// <param name="kernelFP forward projection kernel"></param>
	/// <param name="kernelBP backprojection kernel"></param>
	/// <param name="kernelNLM NLM kernel"></param>
	/// <param name="kernelMed MRP kernel"></param>
	/// <param name="kernelRDP RDP kernel"></param>
	/// <param name="programFP program containing forward projection"></param>
	/// <param name="programBP program containing backprojection"></param>
	/// <param name="programAux program containing NLM/MRP/RDP"></param>
	/// <param name="MethodList reconstruction algorithms selected"></param>
	/// <param name="w_vec specifies some of the special options/parameters used"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters"></param>
	/// <returns></returns>
        inline STATUS_t createKernels(KERNELHANDLE_t & kernelFP, KERNELHANDLE_t & kernelBP, KERNELHANDLE_t & kernelNLM, KERNELHANDLE_t & kernelMed,
			KERNELHANDLE_t & kernelRDP, KERNELHANDLE_t & kernelGGMRF, const PROGRAMHANDLE_t & programFP, const PROGRAMHANDLE_t & programBP, const PROGRAMHANDLE_t & programAux,
#if defined(METAL) || defined(OPENCL)
			const PROGRAMHANDLE_t & programSens, 
#endif // END METAL/OPENCL
            const RecMethods & MethodList, const Weighting & w_vec, const scalarStruct & inputScalars, const int type = -1) {
			STATUS_t status = SUCCESS_VALUE;
#if defined(METAL)
#ifdef AF
			MTL::CommandQueue* arrayFireQueue = afmtl::getQueue();
			queueFP = NS::RetainPtr(arrayFireQueue);
			queueBP = NS::RetainPtr(arrayFireQueue);
#else
			queueFP = NS::TransferPtr(mtlDevice->newCommandQueue());
			queueBP = NS::TransferPtr(mtlDevice->newCommandQueue());
#endif
			if (!queueFP || !queueBP) {
				mexPrint("Unable to create Metal command queues");
				return -1;
			}
#endif // END METAL
			// Kernel for the OS-methods (OSEM, RAMLA, RBI, BSREM, etc.)
			if (inputScalars.FPType == 4 || inputScalars.BPType == 4) {
				if (inputScalars.FPType == 4) {
					CREATE_KERNEL(kernelFP, programFP, "projectorType4Forward", "Failed to create projector type 4 FP kernel\n");
					if (DEBUG || inputScalars.verbose >= 3) {
						mexPrint(BACKEND_STR " kernel for projector type 4 FP successfully created\n");
					}
				}
				if (inputScalars.BPType == 4) {
					if (inputScalars.FPType == 4 && inputScalars.CT)
						GET_KERNEL(kernelBP, programBP, "projectorType4Backward");
					else if (!inputScalars.CT)
						GET_KERNEL(kernelBP, programBP, "projectorType4Forward");
					else
						CREATE_KERNEL(kernelBP, programBP, "projectorType4Backward", "Failed to create projector type 4 BP kernel\n");
					if (DEBUG || inputScalars.verbose >= 3) {
						mexPrint(BACKEND_STR " kernel for projector type 4 BP successfully created\n");
					}
				}
			}
			if (inputScalars.FPType == 5 || inputScalars.BPType == 5) {
				if (inputScalars.FPType == 5) {
					CREATE_KERNEL(kernelFP, programFP, "projectorType5Forward", "Failed to create projector type 5 FP kernel\n");
					if (DEBUG || inputScalars.verbose >= 3) {
						mexPrint(BACKEND_STR " kernel for projector type 5 FP successfully created\n");
					}
				}
				if (inputScalars.BPType == 5) {
					if (inputScalars.FPType == 5)
						GET_KERNEL(kernelBP, programFP, "projectorType5Backward");
					else
						GET_KERNEL(kernelBP, programBP, "projectorType5Backward");
					KCHECK("Failed to create projector type 5 BP kernel\n");
					if (DEBUG || inputScalars.verbose >= 3) {
						mexPrint(BACKEND_STR " kernel for projector type 5 BP successfully created\n");
					}
				}
			}
			if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3 || inputScalars.FPType == 1 || 
				inputScalars.FPType == 2 || inputScalars.FPType == 3) {
				if ((inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3)) {
					GET_KERNEL(kernelFP, programFP, "projectorType123");
					KCHECK("Failed to create projector type 1-3 FP kernel\n");
					if (DEBUG || inputScalars.verbose >= 3) {
						mexPrint(BACKEND_STR " kernel for projector type 1-3 FP successfully created\n");
					}
				}
				if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
					GET_KERNEL(kernelBP, programBP, "projectorType123");
					KCHECK("Failed to create projector type 1-3 BP kernel\n");
					if (DEBUG || inputScalars.verbose >= 3) {
						mexPrint(BACKEND_STR " kernel for projector type 1-3 BP successfully created\n");
					}
				}
			}
#if !defined(METAL)
		if (MethodList.NLM) {
			CREATE_KERNEL(kernelNLM, programAux, "NLM", "Failed to create NLM kernel\n");
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("NLM kernel successfully created\n");
			}
		}
		if (MethodList.MRP) {
			CREATE_KERNEL(kernelMed, programAux, "medianFilter3D", "Failed to create Median kernel\n");
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("Median kernel successfully created\n");
			}
		}
		if (MethodList.RDP) {
			CREATE_KERNEL(kernelRDP, programAux, "RDPKernel", "Failed to create RDP kernel\n");
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("RDP kernel successfully created\n");
			}
		}
		if (MethodList.GGMRF) {
			CREATE_KERNEL(kernelGGMRF, programAux, "GGMRFKernel", "Failed to create GGMRF kernel\n");
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("GGMRF kernel successfully created\n");
			}
		}
		if (MethodList.TV || MethodList.APLS) {
			CREATE_KERNEL(kernelTV, programAux, "TVKernel", "Failed to create TV kernel\n");
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("TV kernel successfully created\n");
			}
		}
		if (MethodList.hyperbolic) {
			CREATE_KERNEL(kernelHyper, programAux, "hyperbolicKernel", "Failed to create hyperbolic prior kernel\n");
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("Hyperbolic prior kernel successfully created\n");
			}
		}
		if (MethodList.PKMA || MethodList.BSREM || MethodList.MBSREM || MethodList.MRAMLA || MethodList.RAMLA) {
			CREATE_KERNEL(kernelPoisson, programAux, "PoissonUpdate", "Failed to create Poisson Update kernel\n");
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("Poisson Update kernel successfully created\n");
			}
		}
#endif // END non-Metal auxiliary kernel creation
		if (MethodList.CPType) {
			CREATE_KERNEL(kernelPDHG, programAux, "PDHGUpdate", "Failed to create PDHG Update kernel\n");
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("PDHG Update kernel successfully created\n");
			}
		}
#if !defined(METAL)
		if (MethodList.ProxTV) {
			GET_KERNEL(kernelProxTVq, programAux, "ProxTVq");
			GET_KERNEL(kernelProxTVDiv, programAux, "ProxTVDivergence");
			GET_KERNEL(kernelProxTVGrad, programAux, "ProxTVGradient");
			CHECK(status, "Failed to create proximal TV kernel\n", status);
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("Proximal TV kernel successfully created\n");
			}
		}
		if (MethodList.ProxRDP) {
			GET_KERNEL(kernelProxq, programAux, "Proxq");
			GET_KERNEL(kernelProxRDP, programAux, "ProxRDP");
			GET_KERNEL(kernelProxTrans, programAux, "ProxTrans");
			KCHECK("Failed to create proximal RDP kernel\n");
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("Proximal RDP kernel successfully created\n");
			}
		}
		if (MethodList.ProxNLM) {
			GET_KERNEL(kernelProxq, programAux, "Proxq");
			GET_KERNEL(kernelProxNLM, programAux, "ProxNLM");
			GET_KERNEL(kernelProxTrans, programAux, "ProxTrans");
			KCHECK("Failed to create proximal NLM kernel\n");
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("Proximal NLM kernel successfully created\n");
			}
		}
		if (MethodList.ProxTGV) {
			GET_KERNEL(kernelProxTVq, programAux, "ProxTVq");
			GET_KERNEL(kernelProxTGVq, programAux, "ProxTGVq");
			GET_KERNEL(kernelProxTVDiv, programAux, "ProxTVDivergence");
			GET_KERNEL(kernelProxTVGrad, programAux, "ProxTVGradient");
			GET_KERNEL(kernelProxTGVDiv, programAux, "ProxTGVDivergence");
			GET_KERNEL(kernelProxTGVSymmDeriv, programAux, "ProxTGVSymmDeriv");
			CHECK(status, "Failed to create proximal TGV kernel\n", status);
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("Proximal TGV kernel successfully created\n");
			}
		}
		if (w_vec.precondTypeMeas[1] || w_vec.precondTypeIm[5]) {
			CREATE_KERNEL(kernelElementMultiply, programAux, "vectorElementMultiply", "Failed to create element-wise multiplication kernel\n");
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("Element-wise kernels successfully created\n");
			}
			CREATE_KERNEL(kernelElementDivision, programAux, "vectorElementDivision", "Failed to create element-wise division kernel\n");
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("Element-wise kernels successfully created\n");
			}
		}
#if defined(OPENCL)
		if (type == 0) {
			kernelsumma = cl::Kernel(programAux, "summa", &status);
			OCL_CHECK(status, "Failed to create implementation 3 kernels\n", -1);
			kernelEstimate = cl::Kernel(programAux, "computeEstimate", &status);
			OCL_CHECK(status, "Failed to create implementation 3 kernels\n", -1);
			kernelForward = cl::Kernel(programAux, "forward", &status);
			if (inputScalars.use_psf) {
				kernelPSFf = cl::Kernel(programAux, "Convolution3D_f", &status);
				kernelPSF = cl::Kernel(programAux, "Convolution3D", &status);
			}
			OCL_CHECK(status, "Failed to create implementation 3 kernels\n", -1);
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("Implementation 3 kernels successfully created\n");
			}
		}
#endif // END CUDA
#endif // END non-Metal auxiliary kernel creation
		if (inputScalars.computeSensImag && inputScalars.listmode > 0) {
			if (inputScalars.BPType == 4)
				GET_KERNEL(kernelSensList, programSens, "projectorType4Forward");
			else
				CREATE_KERNEL(kernelSensList, programSens, "projectorType123", "Failed to create sensitivity image kernels\n");
		}
		if (inputScalars.projector_type == 6) {
			CREATE_KERNEL(kernelRotate, programAux, "rotate", "Failed to create bilinear rotation kernel\n");
			if (DEBUG || inputScalars.verbose >= 3) {
				mexPrint("Bilinear rotation kernel successfully created\n");
			}
		}
		return status;
	}
public:
	ProjectorClass()
#if defined(METAL)
		: autoreleasePool(NS::TransferPtr(NS::AutoreleasePool::alloc()->init()))
#endif
	{}

	inline DEVBUFF_t makeDeviceBuffer(const size_t bytes, const UINT64_t flags, STATUS_t& status) {
		DEVBUFF_t buffer{};
		ALLOC_BUFFER(buffer, flags, bytes);
		return buffer;
	}

	inline STATUS_t writeDeviceBuffer(DEVBUFF_t& buffer, const void* input, const size_t bytes) {
		STATUS_t status = SUCCESS_VALUE;
		WRITE_BUFFER(buffer, bytes, input);
		return status;
	}

	inline STATUS_t readDeviceBuffer(const DEVBUFF_t& buffer, void* output, const size_t bytes) {
		STATUS_t status = SUCCESS_VALUE;
		READ_BUFFER(buffer, bytes, output);
		return status;
	}

	template <typename T>
	inline STATUS_t fillDeviceBuffer(DEVBUFF_t& buffer, const T value, const size_t bytes) {
#if defined(OPENCL)
		return CLCommandQueue[0].enqueueFillBuffer(buffer, value, 0, bytes);
#else
		std::vector<unsigned char> fillData(bytes);
		for (size_t offset = 0; offset < bytes; offset += sizeof(value)) {
			const size_t remaining = bytes - offset;
			const size_t copyBytes = remaining < sizeof(value) ? remaining : sizeof(value);
			std::memcpy(fillData.data() + offset, &value, copyBytes);
		}
		STATUS_t status = SUCCESS_VALUE;
		WRITE_BUFFER(buffer, bytes, fillData.data());
		return status;
#endif
	}

	inline STATUS_t finishDeviceQueue() {
		STATUS_t status = SUCCESS_VALUE;
		FINISH_QUEUE(status, "Queue finish failed\n", status);
		return status;
	}

	// Create and populate a floating-point 3D texture through the backend-specific
	// texture compatibility macro. This is also available to future CUDA callers.
	inline STATUS_t createFloatTexture3DFromHost(TEX3D_t& texture, TEXARRAY_t& array, const float* source,
		const size_t width, const size_t height, const size_t depth) {
		STATUS_t status = SUCCESS_VALUE;
		CREATE_FLOAT_TEXTURE3D_FROM_HOST(texture, array, source, width, height, depth,
			BACKEND_TEXTURE_POINT, BACKEND_TEXTURE_DEFAULT_FLAGS);
		return status;
	}

	inline STATUS_t createFloatTexture3DFromDevice(TEX3D_t& texture, TEXARRAY_t& array, const AFDEVBUFF_t& source,
		const size_t width, const size_t height, const size_t depth) {
		STATUS_t status = SUCCESS_VALUE;
		CREATE_FLOAT_TEXTURE3D_FROM_DEVICE(texture, array, source, width, height, depth,
			BACKEND_TEXTURE_POINT, BACKEND_TEXTURE_DEFAULT_FLAGS);
		return status;
	}

#if defined(METAL)
	NS::SharedPtr<MTL::Device> mtlDevice;
	NS::SharedPtr<MTL::CommandQueue> queueFP, queueBP;
	METAL_im_vectors vec_opencl;
	ScalarKernelParams kParams;
#endif // END METAL
#if defined(CUDA) || defined(HIP)
	std::vector<CUdevice> CUDeviceID;
	std::vector<CUstream> CLCommandQueue;
#if !defined(AF)
	CUcontext standaloneContext = nullptr;
	CUdevice standaloneDevice = 0;
#endif
#elif defined(OPENCL)
	cl::Context CLContext;
	std::vector<cl::Device> CLDeviceID;
	std::vector<cl::CommandQueue> CLCommandQueue;
	OpenCL_im_vectors vec_opencl;
#endif // END CUDA
	KERNELHANDLE_t kernelMBSREM, kernelFP, kernelBP, kernelNLM, kernelMed, kernelRDP, kernelProxTVq, kernelProxTVDiv, kernelProxTVGrad, 
		kernelElementMultiply, kernelElementDivision, kernelTV, kernelProxTGVSymmDeriv, kernelProxTGVDiv, kernelProxTGVq, kernelPoisson, 
		kernelPDHG, kernelProxRDP, kernelProxq, kernelProxTrans, kernelProxNLM, kernelGGMRF, kernelsumma, kernelEstimate, kernelPSF, 
		kernelPSFf, kernelDiv, kernelMult, kernelForward, kernelSensList, kernelApu, kernelHyper, kernelRotate;
	// Device buffers shared across backends
	DEVBUFF_t d_xcenter, d_ycenter, d_zcenter, d_V, d_TOFCenter, d_eFOVIndices, d_weights, d_angle, d_g, d_uref,
		d_rayShiftsDetector, d_rayShiftsSource, d_maskPriorB;
	std::vector<std::vector<DEVBUFF_t>> d_maskBPB;
	std::vector<DEVBUFF_t> d_attenB;
	AFDEVBUFF_t d_output, d_meanBP, d_meanFP, d_inputB, d_W, d_gaussianNLM;
	AFDEVBUFF_t d_qX, d_qY, d_qZ;
	AFDEVBUFF_t d_rX, d_rY, d_rXY, d_rZ, d_rXZ, d_rYZ;
	AFDEVBUFF_t d_vX, d_vY, d_vZ;
	AFDEVBUFF_t d_vector, d_input;
	AFDEVBUFF_t d_im, d_rhs, d_U, d_refIm, d_RDPref;
	// Sensitivity image used by fastPDHG
	AFDEVBUFF_t d_precond;
	AFDEVBUFF_t d_outputCT;
	TEX2D_t d_maskFP, d_maskPrior;
	std::vector<std::vector<TEX2D_t>> d_maskBP;
	std::vector<std::vector<TEX3D_t>> d_maskBP3;
	TEX3D_t d_maskPrior3;
	TEX3D_t d_inputImage, d_urefIm, d_inputI, d_RDPrefI;
	// Texture3D d_imageX, d_imageY; // Unused
	TEXARRAY_t atArray, uRefArray, maskArrayPrior, BPArray, FPArray, integArrayXY, imArray;
	std::vector<std::vector<TEXARRAY_t>> maskArrayBP;
	std::vector<TEX3D_t> d_attenIm;
#if defined(CUDA) || defined(HIP)
	CUmodule programFP, programBP, programAux, programSens;
	std::vector<void*> FPArgs, BPArgs, SensArgs;
	CUDA_im_vectors vec_opencl;
#endif // END CUDA
	std::chrono::steady_clock::time_point tStartLocal, tStartGlobal, tStartAll;
	std::chrono::steady_clock::time_point tEndLocal, tEndGlobal, tEndAll;
	// Distance from the origin to the corner of the image, voxel size and distance from the origin to the opposite corner of the image
	std::vector<FLOAT3_t> b, d, bmax;
	// Axial extent of the full volume with largeDim, where b/bmax above only cover the current subvolume
	float bzGlobalFP[2] = { 0.f, 0.f }, bzGlobalBP[2] = { 0.f, 0.f };
	std::vector<INT3_t> d_N;
	UCHAR_t no_norm = 0;
	size_t memSize = 0ULL;

#if defined(OPENCL)
	// Raw data support back in the future?
//  std::vector<DEVBUFF_t> d_LFull, d_zindexFull, d_xyindexFull; // Unused
	// Image origin
	cl::detail::size_t_array origin = { { 0, 0, 0 } };
	cl::detail::size_t_array region = { { 0, 0, 0 } };
	// Image format
	cl::ImageFormat format;
	cl::ImageFormat formatMask;
#elif defined(METAL)
	std::array<NS::UInteger, 3> origin = { 0, 0, 0 };
	std::array<NS::UInteger, 3> region = { 0, 0, 0 };
#endif // END CUDA
	std::vector<std::vector<TEX3D_t>> d_maskFP3;
	std::vector<std::vector<TEXARRAY_t>> maskArrayFP;
	std::vector<AFDEVBUFF_t> d_Summ;
	std::vector<AFDEVBUFF_t> d_meas, d_rand, d_imTemp, d_imFinal;
	// Vector device buffers common to both backends
	std::vector<std::vector<DEVBUFF_t>> d_maskFPB;
	std::vector<std::vector<DEVBUFF_t>> d_detectorVector;
	std::vector<DEVBUFF_t> d_normFull, d_scatFull, d_xFull, d_zFull;
	std::vector<DEVBUFF_t> d_L;
	std::vector<DEVBUFF_t> d_zindex, d_xyindex, d_T;
	std::vector<std::vector<DEVBUFF_t>> d_norm, d_atten, d_scat, d_x, d_z, d_trIndex, d_axIndex, d_TOFIndex;
	// Precomputed per-projection geometry for the BDD backprojection (16 floats per projection, one buffer per subset)
	std::vector<DEVBUFF_t> d_geomProj5;
	// Host-side storage for the above geometry
	std::vector<std::vector<float>> geomProj5Host;
	std::vector<std::vector<size_t>> erotusBP, erotusPDHG;
#if defined(CUDA) || defined(HIP)
	// This is used to define the additional queues/streams in the multi-resolution case
	// Note that these are only used in the multi-resolution case
	std::vector<CUstream> sideQueues;
	CUevent evMain = nullptr;
	std::vector<CUevent> evSide;
	// Persistent per-volume FP input arrays/textures
	// Previously these were always recreated
	// The downside is that multi-resolution case uses more memory now
	std::vector<CUarray> FPArrayCache;
	std::vector<CUtexObject> FPTexCache;
	// Texture object for fastPDHG when using FPType == 4 --> priors require same texture but integer coordinates
	// with nearest neighbor
	std::vector<CUtexObject> FPTexCachePrior;
	std::vector<CUarray> intArrayCacheXY;
	std::vector<CUtexObject> intTexCacheXY;
	std::vector<CUarray> intArrayCacheYZ;
	std::vector<CUtexObject> intTexCacheYZ;
#elif defined(METAL)
	// One queue and retained last command buffer per multi-resolution volume.
	// Waiting for the last command buffer is sufficient because Metal queues
	// execute their own command buffers in submission order.
	std::vector<NS::SharedPtr<MTL::CommandQueue>> sideQueues;
	std::vector<NS::SharedPtr<MTL::CommandBuffer>> sideCommandBuffers;
	// Persistent per-volume FP input textures, refreshed from the current
	// ArrayFire image estimate before every forward projection.
	std::vector<TEX3D_t> FPTexCache;
#elif defined(OPENCL)
	std::vector<cl::CommandQueue> sideQueues;
	// Persistent per-volume FP input images
	std::vector<cl::Image3D> d_imageCache;
	std::vector<cl::Image3D> d_intImageCacheXY;
	std::vector<cl::Image3D> d_intImageCacheYZ;
#endif
	// Cached dimensions of the above (3 entries per volume; third entry 0 = not yet created)
	std::vector<size_t> imageCacheDims, intImageCacheDimsXY, intImageCacheDimsYZ;
	// Cached dimensions of the persistent BP input image/texture (third entry 0 = not yet created)
	size_t BPImageDims[3] = { 0, 0, 0 };
	// ==== End Mod additions ====
	// fastPDHG requirements
	// Set to 1 only for the backprojections of the actual reconstruction loop
	// The same kernel is also used for the sensitivity image and the power method, where the image estimate must not be touched
	uint8_t fastStep = 0;
	uint8_t fastPositivity = 0;
	// Boolean value to determine whether regularization is used (NLM, RDP, GGMRF, hyperbolic prior, gradient-based TV)
	bool fastNLMUsed = false;
	// 0 = PDHG (and its variants, default), 1 = Poisson (PKMA/MBSREM/BSREM)
	uint8_t fastAlg = 0;
	// Scalar values required specifically for the fastPDHG
	float fastTheta = 0.f, fastTau = 0.f, fastBeta = 0.f, fastEpps = 1e-6f;
	// Poisson-algorithm (PKMA/MBSREM/BSREM) relaxation parameter and alphaM (PKMA) / U (MBSREM) / 1 (BSREM)
	float fastLambda = 0.f, fastAlpha = 1.f;
	// Prior parameters, always passed as a block: h (sigma^2), gamma, GGMRF p/q/c and the adaptive NLM constant
	float fastNLM[6] = { 1.f, 1.f, 0.f, 0.f, 0.f, 0.f };
	// largeDim axial offsets
	// Only for largeDim cases
	int64_t fastImOffset = 0;
	int fastPriorZOffset = 0;
#if defined(CUDA) || defined(HIP)
	~ProjectorClass() {
		if (memAlloc.FPMod)
			getErrorString(cuModuleUnload(programFP));
		if (memAlloc.BPMod)
			getErrorString(cuModuleUnload(programBP));
		if (memAlloc.auxMod)
			getErrorString(cuModuleUnload(programAux));
		if (memAlloc.SensMod)
			getErrorString(cuModuleUnload(programSens));
		if (memAlloc.attenM) {
			for (const auto& timestepBuffers : d_atten)
				for (const auto& buffer : timestepBuffers)
					getErrorString(cuMemFree(buffer));
		}
		if (memAlloc.V)
			getErrorString(cuMemFree(d_V));
		if (memAlloc.atten && !memAlloc.useBuffers) {
			getErrorString(cuArrayDestroy(atArray));
		}
		for (int tt = 0; tt < memAlloc.tSteps; tt++) {
			if (memAlloc.atten && memAlloc.attenSize < tt) {
				if (memAlloc.useBuffers) {
					getErrorString(cuMemFree(d_attenB[tt]));
				}
				else {
					getErrorString(cuTexObjectDestroy(d_attenIm[tt]));
				}
			}
			if (memAlloc.xSteps >= 0) {
				for (int kk = 0; kk <= memAlloc.xSteps / memAlloc.tSteps; kk++) {
					getErrorString(cuMemFree(d_x[tt][kk]));
				}
			}
			if (memAlloc.zType == 0) {
				getErrorString(cuMemFree(d_z[tt][memAlloc.zSteps]));
			}
			else if (memAlloc.zType == 1) {
				for (int kk = 0; kk <= memAlloc.zSteps / memAlloc.tSteps; kk++) {
					getErrorString(cuMemFree(d_z[tt][kk]));
				}
			}
			if (memAlloc.extra) {
				for (int kk = 0; kk < memAlloc.eSteps / memAlloc.tSteps; kk++) {
					getErrorString(cuMemFree(d_scat[tt][kk]));
				}
			}
			if (memAlloc.indexBased) {
				// When the data is loaded one subset at a time, only a single buffer at [0][0] exists
				int nIndex = memAlloc.iSteps / memAlloc.tSteps;
				if (nIndex == 0)
					nIndex = tt == 0 ? memAlloc.iSteps : 0;
				for (int kk = 0; kk < nIndex; kk++) {
					getErrorString(cuMemFree(d_trIndex[tt][kk]));
					getErrorString(cuMemFree(d_axIndex[tt][kk]));
				}
			}
			if (memAlloc.TOFIndex) {
				// When the data is loaded one subset at a time, only a single buffer at [0][0] exists
				int nTOF = memAlloc.TOFSteps / memAlloc.tSteps;
				if (nTOF == 0)
					nTOF = tt == 0 ? memAlloc.TOFSteps : 0;
				for (int kk = 0; kk < nTOF; kk++) {
					getErrorString(cuMemFree(d_TOFIndex[tt][kk]));
				}
			}
		}
		if (memAlloc.offsetT) {
			for (int kk = 0; kk < memAlloc.oSteps; kk++) {
				getErrorString(cuMemFree(d_T[kk]));
			}
		}
		if (memAlloc.geom5) {
			for (int kk = 0; kk < memAlloc.g5Steps; kk++) {
				getErrorString(cuMemFree(d_geomProj5[kk]));
			}
		}
		if (memAlloc.TOF) {
			getErrorString(cuMemFree(d_TOFCenter));
		}
		if (memAlloc.eFOV) {
			getErrorString(cuMemFree(d_eFOVIndices));
		}
		if (memAlloc.rayShifts) {
			getErrorString(cuMemFree(d_rayShiftsDetector));
			getErrorString(cuMemFree(d_rayShiftsSource));
			for (auto& perTimestep : d_detectorVector)
				for (auto& buffer : perTimestep)
					getErrorString(cuMemFree(buffer));
		}
		if (memAlloc.GGMRF) {
			getErrorString(cuMemFree(d_weights));
		}
		if (memAlloc.norm) {
			for (const auto& timestepBuffers : d_norm)
				for (const auto& buffer : timestepBuffers)
					getErrorString(cuMemFree(buffer));
		}
		if (memAlloc.angle)
			getErrorString(cuMemFree(d_angle));
		if (memAlloc.xFull)
			getErrorString(cuMemFree(d_xFull[0]));
		if (memAlloc.zFull)
			getErrorString(cuMemFree(d_zFull[0]));
		if (memAlloc.maskFP) {
			if (memAlloc.useBuffers) {
				for (const auto& timestepBuffers : d_maskFPB)
					for (const auto& buffer : timestepBuffers)
						getErrorString(cuMemFree(buffer));
			}
			else {
				if (d_maskFP3.size() > 0) {
					for (const auto& timestepTextures : d_maskFP3)
						for (const auto& texture : timestepTextures)
							getErrorString(cuTexObjectDestroy(texture));
				}
				else {
					getErrorString(cuTexObjectDestroy(d_maskFP));
				}
				for (const auto& timestepArrays : maskArrayFP)
					for (const auto& array : timestepArrays)
						getErrorString(cuArrayDestroy(array));
			}
		}
		if (memAlloc.maskBP) {
			if (memAlloc.useBuffers) {
				for (const auto& timestepBuffers : d_maskBPB)
					for (const auto& buffer : timestepBuffers)
						getErrorString(cuMemFree(buffer));
			}
			else {
				if (!d_maskBP3.empty()) {
					for (const auto& timestepTextures : d_maskBP3)
						for (const auto& texture : timestepTextures)
							getErrorString(cuTexObjectDestroy(texture));
				}
				else {
					for (const auto& timestepTextures : d_maskBP)
						for (const auto& texture : timestepTextures)
							getErrorString(cuTexObjectDestroy(texture));
				}
				for (const auto& timestepArrays : maskArrayBP)
					for (const auto& array : timestepArrays)
						getErrorString(cuArrayDestroy(array));
			}
		}
		if (memAlloc.priorMask) {
			if (memAlloc.useBuffers) {
				getErrorString(cuMemFree(d_maskPriorB));
			}
			else {
				getErrorString(cuTexObjectDestroy(d_maskPrior));
				getErrorString(cuArrayDestroy(maskArrayPrior));
			}
		}
		if (memAlloc.NLMRef == 1) {
			getErrorString(cuTexObjectDestroy(d_urefIm));
			getErrorString(cuArrayDestroy(uRefArray));
		}
		else if (memAlloc.NLMRef == 2) {
			getErrorString(cuMemFree(d_uref));
		}
		for (size_t kk = 0; kk < FPTexCache.size(); kk++) {
			if (imageCacheDims[kk * 3 + 2] != 0) {
				getErrorString(cuTexObjectDestroy(FPTexCache[kk]));
				if (kk < FPTexCachePrior.size())
					getErrorString(cuTexObjectDestroy(FPTexCachePrior[kk]));
				getErrorString(cuArrayDestroy(FPArrayCache[kk]));
			}
		}
		for (size_t kk = 0; kk < intTexCacheXY.size(); kk++) {
			if (intImageCacheDimsXY[kk * 3 + 2] != 0) {
				getErrorString(cuTexObjectDestroy(intTexCacheXY[kk]));
				getErrorString(cuArrayDestroy(intArrayCacheXY[kk]));
			}
		}
		for (size_t kk = 0; kk < intTexCacheYZ.size(); kk++) {
			if (intImageCacheDimsYZ[kk * 3 + 2] != 0) {
				getErrorString(cuTexObjectDestroy(intTexCacheYZ[kk]));
				getErrorString(cuArrayDestroy(intArrayCacheYZ[kk]));
			}
		}
		if (BPImageDims[2] != 0) {
			getErrorString(cuTexObjectDestroy(d_inputImage));
			getErrorString(cuArrayDestroy(BPArray));
		}
		for (size_t kk = 0; kk < sideQueues.size(); kk++)
			getErrorString(cuStreamDestroy(sideQueues[kk]));
		for (size_t kk = 0; kk < evSide.size(); kk++)
			getErrorString(cuEventDestroy(evSide[kk]));
		if (evMain != nullptr)
			getErrorString(cuEventDestroy(evMain));
#if !defined(AF)
		if (!CLCommandQueue.empty())
			getErrorString(cuStreamDestroy(CLCommandQueue[0]));
		if (standaloneContext != nullptr)
			getErrorString(cuDevicePrimaryCtxRelease(standaloneDevice));
#endif
	}
#elif defined(OPENCL) || defined(METAL)
	~ProjectorClass() {}
#endif // END CUDA

	/// <summary>
	/// Create the additional queues/streams for the multi-resolution case,
	//// The main queue/stream is always the ArrayFire queue/stream, the other volumes use their own ones
	/// </summary>
	/// <param name="n">Number of queues, one for the main image plus each multi-resolution volume</param>
	inline int initSideQueues(const int n) {
		if (n <= 0 || static_cast<int>(sideQueues.size()) >= n)
			return 0;
#if defined(METAL)
		for (int kk = static_cast<int>(sideQueues.size()); kk < n; kk++) {
			auto queue = NS::TransferPtr(mtlDevice->newCommandQueue());
			if (!queue) {
				mexPrint("Failed to create Metal side command queue\n");
				return -1;
			}
			sideQueues.push_back(std::move(queue));
			sideCommandBuffers.emplace_back(nullptr);
		}
		return 0;
#else
#if defined(CUDA) || defined(HIP)
		CUresult status = CUDA_SUCCESS;
		if (evMain == nullptr) {
			status = cuEventCreate(&evMain, CU_EVENT_DISABLE_TIMING);
			CUDA_CHECK(status, "Failed to create main stream event\n", -1);
		}
		for (int kk = static_cast<int>(sideQueues.size()); kk < n; kk++) {
			CUstream s;
			// Create the stream
			status = cuStreamCreate(&s, CU_STREAM_NON_BLOCKING);
			CUDA_CHECK(status, "Failed to create side stream\n", -1);
			sideQueues.push_back(s);
			CUevent e;
			// Create the event to make sure all relevant computations are complete after the BP
			status = cuEventCreate(&e, CU_EVENT_DISABLE_TIMING);
			CUDA_CHECK(status, "Failed to create side stream event\n", -1);
			evSide.push_back(e);
		}
#else
		cl_int status = CL_SUCCESS;
		for (int kk = static_cast<int>(sideQueues.size()); kk < n; kk++) {
			sideQueues.push_back(cl::CommandQueue(CLContext, CLDeviceID[0], 0, &status));
			OCL_CHECK(status, "Failed to create side command queue\n", -1);
		}
#endif
		return 0;
#endif
	}

	/// <summary>
	/// This function forces the main queue/stream to wait for all the side queues/streams
	// This makes sure that the computations will only proceed after all the BPs are done
	/// </summary>
	inline int joinSideQueues() {
#if defined(CUDA) || defined(HIP)
		CUresult status = CUDA_SUCCESS;
		for (size_t kk = 0; kk < sideQueues.size(); kk++) {
			status = cuEventRecord(evSide[kk], sideQueues[kk]);
			CUDA_CHECK(status, "Failed to record side stream event\n", -1);
			status = cuStreamWaitEvent(CLCommandQueue[0], evSide[kk], 0);
			CUDA_CHECK(status, "Failed to make main stream wait for side stream\n", -1);
		}
#elif defined(METAL)
		int result = 0;
		for (auto& commandBuffer : sideCommandBuffers) {
			if (!commandBuffer)
				continue;
			commandBuffer->waitUntilCompleted();
			if (commandBuffer->status() == MTL::CommandBufferStatusError) {
				NS::Error* error = commandBuffer->error();
				const char* message = error && error->localizedDescription()
					? error->localizedDescription()->utf8String() : "unknown Metal error";
				mexPrintBase("Metal side-queue backprojection failed: %s\n", message);
				result = -1;
			}
			commandBuffer.reset();
		}
		return result;
#elif defined(OPENCL)
		cl_int status = CL_SUCCESS;
		if (sideQueues.size() > 0) {
			std::vector<cl::Event> waitList(sideQueues.size());
			for (size_t kk = 0; kk < sideQueues.size(); kk++) {
				status = sideQueues[kk].enqueueMarkerWithWaitList(nullptr, &waitList[kk]);
				OCL_CHECK(status, "Failed to enqueue side queue marker\n", -1);
				status = sideQueues[kk].flush();
				OCL_CHECK(status, "Failed to flush side queue\n", -1);
			}
			status = CLCommandQueue[0].enqueueBarrierWithWaitList(&waitList);
			OCL_CHECK(status, "Failed to enqueue main queue barrier\n", -1);
		}
#endif
		return 0;
	}

#ifndef METAL
#if defined(CUDA) || defined(HIP)
	// Emits the build options for fastPDHG. nVoxZ is the number of voxels one work-item handles in z
	inline void addFastPDHGOptions(std::vector<std::string> & os_options, const scalarStruct & inputScalars, const Weighting & w_vec,
		const RecMethods & MethodList, const size_t nVoxZ) {
#else
	inline void addFastPDHGOptions(std::string & os_options, const scalarStruct & inputScalars, const Weighting & w_vec,
		const RecMethods & MethodList, const size_t nVoxZ) {
#endif
		ADD_OPT(os_options, "-DFASTPDHG");
		// Determine the algorithm and set the proper definitions
		const bool poissonAlg = MethodList.PKMA || MethodList.MBSREM || MethodList.BSREM;
		if (poissonAlg) {
			fastAlg = 1;
			ADD_OPT_INT(os_options, "-DFASTALG", 1);
			if (MethodList.PKMA) {
				ADD_OPT(os_options, "-DPKMA");
			}
			else if (MethodList.MBSREM) {
				ADD_OPT(os_options, "-DMBSREM");
			}
			else if (MethodList.BSREM) {
				ADD_OPT(os_options, "-DBSREM");
			}
		}
		else
			fastAlg = 0;
		// BSREM does preconditioning only after all subiterations are complete
		if (!MethodList.BSREM) {
		if (w_vec.precondTypeIm[0]) {
			ADD_OPT_INT(os_options, "-DFASTPRECOND", 1);
		}
		else if (w_vec.precondTypeIm[1]) {
			ADD_OPT_INT(os_options, "-DFASTPRECOND", 2);
		}
		// Whether the NLM exclusive local memory caching with patch window is used or not
		bool nlmTileZEmitted = false;
		// TODO: Decrease NVOXELS if the cache is too big?
		auto addFastSWLocalCache = [&](const bool anatomical) {
			if (!FASTNLMLOCALCACHE || anatomical)
				return;
			const size_t cacheSize = (local_size[0] + 2ULL * w_vec.Ndx) * (local_size[1] + 2ULL * w_vec.Ndy)
				* (nVoxZ + 2ULL * w_vec.Ndz) * sizeof(float);
			ADD_OPT(os_options, "-DFASTSWLOCAL");
			if (!nlmTileZEmitted) {
				ADD_OPT_INT(os_options, "-DFASTNLMTILEZ", nVoxZ);
				nlmTileZEmitted = true;
			}
			if (DEBUG) {
				mexPrintBase("Caching the search-window neighborhood of fastPDHG in local memory, %u bytes per work-group\n", static_cast<uint32_t>(cacheSize));
				mexEval();
			}
		};
		// The NLM prior is only included when it was actually selected
		if (MethodList.NLM) {
			fastNLMUsed = true;
			ADD_OPT(os_options, "-DFASTNLM");
			if (w_vec.NLM_MRP) {
				ADD_OPT_INT(os_options, "-DNLTYPE", 2);
			}
			else if (w_vec.NLTV) {
				ADD_OPT_INT(os_options, "-DNLTYPE", 1);
			}
			else if (w_vec.NLRD) {
				ADD_OPT_INT(os_options, "-DNLTYPE", 3);
			}
			else if (w_vec.NLLange) {
				ADD_OPT_INT(os_options, "-DNLTYPE", 4);
			}
			else if (w_vec.NLLangeFiltered) {
				ADD_OPT_INT(os_options, "-DNLTYPE", 5);
			}
			else if (w_vec.NLGGMRF) {
				ADD_OPT_INT(os_options, "-DNLTYPE", 6);
			}
			else {
				ADD_OPT_INT(os_options, "-DNLTYPE", 0);
			}
			if (w_vec.NLAdaptive)
				ADD_OPT(os_options, "-DNLMADAPTIVE");
			if (w_vec.NLM_anatomical)
				ADD_OPT(os_options, "-DNLMREF");
			ADD_OPT_INT(os_options, "-DSWINDOWX", w_vec.Ndx);
			ADD_OPT_INT(os_options, "-DSWINDOWY", w_vec.Ndy);
			ADD_OPT_INT(os_options, "-DSWINDOWZ", w_vec.Ndz);
			ADD_OPT_INT(os_options, "-DPWINDOWX", w_vec.Nlx);
			ADD_OPT_INT(os_options, "-DPWINDOWY", w_vec.Nly);
			ADD_OPT_INT(os_options, "-DPWINDOWZ", w_vec.Nlz);
			// TODO: Decrease NVOXELS if the cache is too big?
			if (FASTNLMLOCALCACHE) {
				const size_t cacheSize = (local_size[0] + 2ULL * (w_vec.Ndx + w_vec.Nlx))
					* (local_size[1] + 2ULL * (w_vec.Ndy + w_vec.Nly))
					* (nVoxZ + 2ULL * (w_vec.Ndz + w_vec.Nlz)) * sizeof(float) * (w_vec.NLM_anatomical ? 2ULL : 1ULL);
				ADD_OPT(os_options, "-DFASTNLMLOCAL");
				ADD_OPT_INT(os_options, "-DFASTNLMTILEZ", nVoxZ);
				nlmTileZEmitted = true;
				if (DEBUG) {
					mexPrintBase("Caching the NLM neighborhood of fastPDHG in local memory, %u bytes per work-group\n", static_cast<uint32_t>(cacheSize));
					mexEval();
				}
			}
		}
		// The RDP prior is only included when it was actually selected
		// RDP reads its neighborhood from the FP cache texture, or, for the "with corners" variant, 
		// optionally from a local-memory cache.
		else if (MethodList.RDP) {
			fastNLMUsed = true;
			ADD_OPT(os_options, "-DFASTRDP");
			if (w_vec.RDPLargeNeighbor) {
				ADD_OPT(os_options, "-DFASTRDPCORNERS");
				ADD_OPT_INT(os_options, "-DSWPRIORTYPE", 0);
				ADD_OPT_INT(os_options, "-DSWINDOWX", w_vec.Ndx);
				ADD_OPT_INT(os_options, "-DSWINDOWY", w_vec.Ndy);
				ADD_OPT_INT(os_options, "-DSWINDOWZ", w_vec.Ndz);
				addFastSWLocalCache(w_vec.RDP_anatomical);
			}
			if (w_vec.RDP_anatomical)
				ADD_OPT(os_options, "-DPRIORREF");
		}
		// Same as above with corners
		else if (MethodList.GGMRF) {
			fastNLMUsed = true;
			ADD_OPT(os_options, "-DFASTGGMRF");
			ADD_OPT_INT(os_options, "-DSWPRIORTYPE", 1);
			ADD_OPT_INT(os_options, "-DSWINDOWX", w_vec.Ndx);
			ADD_OPT_INT(os_options, "-DSWINDOWY", w_vec.Ndy);
			ADD_OPT_INT(os_options, "-DSWINDOWZ", w_vec.Ndz);
			addFastSWLocalCache(false);
		}
		// Same as above
		else if (MethodList.hyperbolic) {
			fastNLMUsed = true;
			ADD_OPT(os_options, "-DFASTHYPER");
			ADD_OPT_INT(os_options, "-DSWPRIORTYPE", 2);
			ADD_OPT_INT(os_options, "-DSWINDOWX", w_vec.Ndx);
			ADD_OPT_INT(os_options, "-DSWINDOWY", w_vec.Ndy);
			ADD_OPT_INT(os_options, "-DSWINDOWZ", w_vec.Ndz);
			addFastSWLocalCache(false);
		}
		// TV includes a lot of different variants, but only the ones without reference images
		// are supported by fastPDHG
		// These include "normal" TV type and TV types 2 and 4
		else if (MethodList.TV) {
			fastNLMUsed = true;
			ADD_OPT(os_options, "-DFASTTV");
			if (w_vec.data.TVtype == 4)
				ADD_OPT(os_options, "-DSATV");
			else if (w_vec.data.TVtype == 2)
				ADD_OPT(os_options, "-DJPTV");
		}
		} // end !MethodList.BSREM
	}
#endif // END METAL

	/// <summary>
	/// This function creates the projector class object
	/// </summary>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="w_vec specifies some of the special options/parameters used"></param>
	/// <param name="MethodList specifies the algorithms and priors used"></param>
	/// <param name="header_directory the location of the kernel and header files"></param>
	/// <returns></returns>
	inline int addProjector(scalarStruct & inputScalars, Weighting & w_vec, const RecMethods & MethodList, const char* header_directory, const int type = -1) {
		// Set-up the local group size
#if defined(CUDA)
		local_size[0] = 32ULL;
#elif defined(METAL) || defined(OPENCL) || defined(HIP)
		// HIP needs to use 64 since CDNA is wave64, not 32 like RDNA
		local_size[0] = 64ULL;
#endif // END CUDA
		local_size[1] = 1ULL;
		local_size[2] = 1ULL;
		if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3 || (inputScalars.BPType == 4 && 
			(!inputScalars.CT && !inputScalars.PET && !inputScalars.SPECT)))
			local_size[0] = 128ULL;
		if (inputScalars.BPType == 4 || inputScalars.BPType == 5 || ((inputScalars.PET || inputScalars.SPECT || inputScalars.CT) && inputScalars.listmode == 0)) {
			if (inputScalars.nColsD > 1 && !(inputScalars.BPType == 4 && (!inputScalars.CT && !inputScalars.PET && !inputScalars.SPECT))) {
				local_size[0] = 16ULL;
				local_size[1] = 16ULL;
			}
		}
		// Override the above defaults with user-supplied values (per dimension). A negative input value
		// keeps the corresponding default computed above.
		for (int ls = 0; ls < 3; ls++)
			if (inputScalars.localSize[ls] > 0)
				local_size[ls] = static_cast<size_t>(inputScalars.localSize[ls]);
		if (DEBUG) {
			mexPrintBase("inputScalars.nColsD = %u\n", inputScalars.nColsD);
			mexPrintBase("inputScalars.nRowsD = %u\n", inputScalars.nRowsD);
			mexPrintBase("local_size[0] = %u\n", local_size[0]);
			mexPrintBase("local_size[1] = %u\n", local_size[1]);
			mexEval();
		}
		// Local group for priors
		local_sizePrior[0] = 16ULL;
		local_sizePrior[1] = 16ULL;
		local_sizePrior[2] = 1ULL;
		STATUS_t status = SUCCESS_VALUE;
#if defined(CUDA) || defined(HIP)
		nvrtcResult status2 = NVRTC_SUCCESS;
#endif // END CUDA

#if defined(CUDA) || defined(HIP)
		// Create the CUDA/HIP context and stream and assign the device
#if defined(AF)
		int af_id = af::getDevice();
		CUDeviceID.push_back(afcu::getNativeId(af_id));
		CLCommandQueue.push_back(afcu::getStream(CUDeviceID[0]));
#else
		status = cuInit(0);
		CUDA_CHECK(status, "Failed to initialize CUDA\n", -1);
		status = cuDeviceGet(&standaloneDevice, inputScalars.platform);
		CUDA_CHECK(status, "Failed to select the CUDA device\n", -1);
		status = cuDevicePrimaryCtxRetain(&standaloneContext, standaloneDevice);
		CUDA_CHECK(status, "Failed to retain the CUDA context\n", -1);
		status = cuCtxSetCurrent(standaloneContext);
		CUDA_CHECK(status, "Failed to activate the CUDA context\n", -1);
		CUDeviceID.push_back(standaloneDevice);
		CUstream stream = nullptr;
		status = cuStreamCreate(&stream, CU_STREAM_DEFAULT);
		CUDA_CHECK(status, "Failed to create the CUDA stream\n", -1);
		CLCommandQueue.push_back(stream);
#endif

		status2 = createProgram(programFP, programBP, programAux, header_directory, inputScalars, MethodList, w_vec, local_size, type);
		if (status2 != NVRTC_SUCCESS) {
			std::cerr << "Error while creating program" << std::endl;
			return -1;
		}
#elif defined(METAL)
		PROGRAMHANDLE_t programFP, programBP, programAux, programSens;
		status = createProgram(programFP, programBP, programAux, programSens, header_directory, inputScalars, MethodList, w_vec, local_size, type);
		CHECK(status, "Error while creating Metal program\n", -1);
#elif defined(OPENCL)
		// Create the OpenCL context and command queue and assign the device
#ifdef AF
		CLContext = afcl::getContext(true);
		CLCommandQueue.push_back(cl::CommandQueue(afcl::getQueue(true), true));
		// Use AF command queue to query for the current OpenCL device
		CLDeviceID.push_back(CLCommandQueue[0].getInfo<CL_QUEUE_DEVICE>(&status));
		OCL_CHECK(status, "\n", -1);
#else
		status = clGetPlatformsContext(inputScalars.platform, CLContext, CLCommandQueue, inputScalars.usedDevices, CLDeviceID);
#endif
		// Force OpenCL CPUs to use buffers
		// Images tend to be much slower on CPUs
		const cl_device_type devType = CLDeviceID[0].getInfo<CL_DEVICE_TYPE>(&status);
		if ((devType & CL_DEVICE_TYPE_CPU) && inputScalars.useImages) {
			inputScalars.useImages = false;
		}
		// For NVIDIA cards, 32 local size seems more optimal with 1D kernelFP (unless the user gave an explicit value)
		std::string deviceName = CLDeviceID[0].getInfo<CL_DEVICE_VENDOR>(&status);
		std::string NV("NVIDIA Corporation");
		if (inputScalars.localSize[0] <= 0 && NV.compare(deviceName) == 0 && (inputScalars.projector_type == 1 || inputScalars.projector_type == 11) 
			&& local_size[1] == 1ULL)
			local_size[0] = 32ULL;
		if (DEBUG) {
			std::string deviceName2 = CLDeviceID[0].getInfo<CL_DEVICE_NAME>(&status);
			UINT64_t apu = CLDeviceID[0].getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>(&status);
			UINT32_t apu2 = CLDeviceID[0].getInfo<CL_DEVICE_ADDRESS_BITS>(&status);
			mexPrintBase("CL_DEVICE_MAX_MEM_ALLOC_SIZE = %llu\n", apu);
			mexPrintBase("CL_DEVICE_ADDRESS_BITS = %u\n", apu2);
			mexPrint(deviceName.c_str());
			mexPrint(deviceName2.c_str());
			mexEval();
		}
#endif // END CUDA
		// Use the prior local sizes
		// TODO: Make the prior local size adjustable
#ifndef METAL
		if (inputScalars.fastPDHG && (inputScalars.BPType == 4 || inputScalars.BPType == 5) && inputScalars.CT && inputScalars.listmode == 0) {
			local_size[0] = local_sizePrior[0];
			local_size[1] = local_sizePrior[1];
			local_size[2] = 1ULL;
		}
		// Check whether the prior local neighborhood can be cached in local memory
		// If not, disable fastPDHG
		if (inputScalars.fastPDHG && FASTNLMLOCALCACHE && inputScalars.CT && inputScalars.listmode == 0 && (inputScalars.BPType == 4 || inputScalars.BPType == 5)) {
			const size_t nVoxZ = (inputScalars.BPType == 5) ? (inputScalars.pitch ? 1ULL : static_cast<size_t>(NVOXELS5))
				: (inputScalars.useHelical ? static_cast<size_t>(NVOXELSHELICAL) : static_cast<size_t>(NVOXELS));
			size_t cacheSize = 0ULL;
			int memLoc = 0;
#if defined(CUDA) || defined(HIP)
			cuDeviceGetAttribute(&memLoc, CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK, CUDeviceID[0]);
#elif defined(OPENCL)
			memLoc = static_cast<int>(CLDeviceID[0].getInfo<CL_DEVICE_LOCAL_MEM_SIZE>(&status));
#endif // END CUDA
			if (MethodList.NLM)
				cacheSize = (local_size[0] + 2ULL * (w_vec.Ndx + w_vec.Nlx)) * (local_size[1] + 2ULL * (w_vec.Ndy + w_vec.Nly))
					* (nVoxZ + 2ULL * (w_vec.Ndz + w_vec.Nlz)) * sizeof(float) * (w_vec.NLM_anatomical ? 2ULL : 1ULL);
			else if ((MethodList.RDP && w_vec.RDPLargeNeighbor && !w_vec.RDP_anatomical) || MethodList.GGMRF || MethodList.hyperbolic)
				cacheSize = (local_size[0] + 2ULL * w_vec.Ndx) * (local_size[1] + 2ULL * w_vec.Ndy)
					* (nVoxZ + 2ULL * w_vec.Ndz) * sizeof(float);
			if (cacheSize >= memLoc - 1024ULL) {
				inputScalars.fastPDHG = false;
				if (DEBUG) {
					mexPrintBase("fastPDHG disabled: the fused prior's local-memory tile needs %u bytes per work-group, which exceeds the limit. Using the regular (multi-step) path.\n", static_cast<uint32_t>(cacheSize));
					mexEval();
				}
			}
		}
#endif
#if defined(CUDA) || defined(HIP)
		if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("CUDA programs successfully created\n");
		}
#elif defined(METAL)
		if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Metal programs successfully created\n");
		}
#elif defined(OPENCL)
		UINT64_t constantBufferSize = CLDeviceID[0].getInfo<CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE>(&status);

		if ((inputScalars.size_of_x + inputScalars.size_z) * sizeof(float) >= constantBufferSize)
			constantBuffer = true;
		if (DEBUG) {
			mexPrintBase("CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE = %u\n", constantBufferSize);
			mexPrintBase("(inputScalars.size_of_x + inputScalars.size_z) * sizeof(float) = %u\n", (inputScalars.size_of_x + inputScalars.size_z) * sizeof(float));
			mexPrintBase("inputScalars.size_of_x = %u\n", inputScalars.size_of_x);
			mexPrintBase("inputScalars.size_z = %u\n", inputScalars.size_z);
			mexEval();
		}
#endif // END CUDA
        if (DEBUG || inputScalars.verbose >= 3)
            mexPrint(BACKEND_STR " programs successfully created\n");

#if defined(CUDA) || defined(HIP)
		status = createKernels(kernelFP, kernelBP, kernelNLM, kernelMed, kernelRDP, kernelGGMRF, programFP, programBP, programAux, 
			MethodList, w_vec, inputScalars, type);
		CUDA_CHECK(status, "Failed to create kernels\n", -1);
#elif defined(METAL)
		status = createKernels(kernelFP, kernelBP, kernelNLM, kernelMed, kernelRDP, kernelGGMRF, programFP, programBP, programAux, 
			programSens, MethodList, w_vec, inputScalars, type);
		CHECK(status, "Failed to create Metal kernels\n", -1);
#elif defined(OPENCL)
		cl::Program programFP, programBP, programAux, programSens;

		status = createProgram(CLContext, CLDeviceID[0], programFP, programBP, programAux, programSens, header_directory, inputScalars, 
			MethodList, w_vec, local_size, type);
		OCL_CHECK(status, "Error while creating program\n", -1);
#endif // END CUDA
		if (DEBUG || inputScalars.verbose >= 3)
			mexPrint(BACKEND_STR " kernels successfully created\n");
#if defined(OPENCL)
		status = createKernels(kernelFP, kernelBP, kernelNLM, kernelMed, kernelRDP, kernelGGMRF, programFP, programBP, programAux, 
			programSens, MethodList, w_vec, inputScalars, type);
		OCL_CHECK(status, "Failed to create kernels\n", -1);
		if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("OpenCL kernels successfully created\n");
		}
		format.image_channel_order = CL_A;
		format.image_channel_data_type = CL_FLOAT;
		formatMask.image_channel_order = CL_A;
		formatMask.image_channel_data_type = CL_UNSIGNED_INT8;
#endif // END CUDA

		if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0) {
			erotus[0] = inputScalars.nRowsD % local_size[0];
			if (inputScalars.FPType == 5)
				erotus[1] = ((inputScalars.nColsD + NVOXELSFP - 1) / NVOXELSFP) % local_size[1];
			else
				erotus[1] = inputScalars.nColsD % local_size[1];
			if (erotus[1] > 0)
				erotus[1] = (local_size[1] - erotus[1]);
			if (erotus[0] > 0)
				erotus[0] = (local_size[0] - erotus[0]);
		}

		if ((MethodList.ProxTGV || MethodList.ProxTV || MethodList.ProxRDP)) {
			erotusPriorEFOV[0] = inputScalars.NxPrior % local_sizePrior[0];
			erotusPriorEFOV[1] = inputScalars.NyPrior % local_sizePrior[1];
			erotusPriorEFOV[2] = inputScalars.NzPrior % local_sizePrior[2];
			if (erotusPriorEFOV[0] > 0)
				erotusPriorEFOV[0] = (local_sizePrior[0] - erotusPriorEFOV[0]);
			if (erotusPriorEFOV[1] > 0)
				erotusPriorEFOV[1] = (local_sizePrior[1] - erotusPriorEFOV[1]);
			if (erotusPriorEFOV[2] > 0)
				erotusPriorEFOV[2] = (local_sizePrior[1] - erotusPriorEFOV[2]);
			SET_LAUNCH_RANGE3(globalPriorEFOV,
				inputScalars.NxPrior + erotusPriorEFOV[0],
				inputScalars.NyPrior + erotusPriorEFOV[1],
				inputScalars.NzPrior + erotusPriorEFOV[2],
				local_sizePrior);
		}

		erotusBP.resize(2);
		erotusPDHG.resize(2);
		if (MethodList.CPType || MethodList.PKMA || MethodList.BSREM || MethodList.RAMLA || MethodList.MRAMLA || MethodList.MBSREM) {
			for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
				erotusPDHG[0].emplace_back(inputScalars.Nx[ii] % local_sizePrior[0]);
				erotusPDHG[1].emplace_back(inputScalars.Ny[ii] % local_sizePrior[1]);
				if (erotusPDHG[0][ii] > 0)
					erotusPDHG[0][ii] = (local_sizePrior[0] - erotusPDHG[0][ii]);
				if (erotusPDHG[1][ii] > 0)
					erotusPDHG[1][ii] = (local_sizePrior[1] - erotusPDHG[1][ii]);
			}
		}
		for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
			erotusBP[0].emplace_back(inputScalars.Nx[ii] % local_size[0]);
			erotusBP[1].emplace_back(inputScalars.Ny[ii] % local_size[1]);
			if (erotusBP[0][ii] > 0)
				erotusBP[0][ii] = (local_size[0] - erotusBP[0][ii]);
			if (erotusBP[1][ii] > 0)
				erotusBP[1][ii] = (local_size[1] - erotusBP[1][ii]);
		}
		SET_RANGE2(local, local_size[0], local_size[1]);
		SET_RANGE3(localPrior, local_sizePrior[0], local_sizePrior[1], local_sizePrior[2]);
		erotusPrior[0] = inputScalars.Nx[0] % local_sizePrior[0];
		erotusPrior[1] = inputScalars.Ny[0] % local_sizePrior[1];
		erotusPrior[2] = inputScalars.Nz[0] % local_sizePrior[2];
		if (erotusPrior[0] > 0)
			erotusPrior[0] = (local_sizePrior[0] - erotusPrior[0]);
		if (erotusPrior[1] > 0)
			erotusPrior[1] = (local_sizePrior[1] - erotusPrior[1]);
		if (erotusPrior[2] > 0)
			erotusPrior[2] = (local_sizePrior[1] - erotusPrior[2]);
		SET_LAUNCH_RANGE3(globalPrior,
			inputScalars.Nx[0] + erotusPrior[0],
			inputScalars.Ny[0] + erotusPrior[1],
			inputScalars.Nz[0] + erotusPrior[2],
			localPrior);
		d_NOrig = make_vec3<INT3_t>(static_cast<INT32_t>(inputScalars.NxOrig), static_cast<INT32_t>(inputScalars.NyOrig), static_cast<INT32_t>(inputScalars.NzOrig));
		d_NPrior = make_vec3<INT3_t>(static_cast<INT32_t>(inputScalars.NxPrior), static_cast<INT32_t>(inputScalars.NyPrior), static_cast<INT32_t>(inputScalars.NzPrior));
		dPitch = { w_vec.dPitchX, w_vec.dPitchY };
		if (inputScalars.SPECT) {
			ellipseCenter = make_vec3<FLOAT3_t>(inputScalars.ellipseCenterX, inputScalars.ellipseCenterY, inputScalars.ellipseCenterZ);
			ellipseRadii = make_vec3<FLOAT3_t>(inputScalars.ellipseRadiusX, inputScalars.ellipseRadiusY, inputScalars.ellipseRadiusZ);
		}
		b.resize(inputScalars.nMultiVolumes + 1);
		d.resize(inputScalars.nMultiVolumes + 1);
		d_N.resize(inputScalars.nMultiVolumes + 1);
		bmax.resize(inputScalars.nMultiVolumes + 1);
		for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
			b[ii] = make_vec3<FLOAT3_t>(inputScalars.bx[ii], inputScalars.by[ii], inputScalars.bz[ii]);
			d[ii] = make_vec3<FLOAT3_t>(inputScalars.dx[ii], inputScalars.dy[ii], inputScalars.dz[ii]);
			d_N[ii] = make_vec3<INT3_t>(static_cast<INT32_t>(inputScalars.Nx[ii]), static_cast<INT32_t>(inputScalars.Ny[ii]), static_cast<INT32_t>(inputScalars.Nz[ii]));
			bmax[ii] = make_vec3<FLOAT3_t>(static_cast<float>(inputScalars.Nx[ii]) * inputScalars.dx[ii] + inputScalars.bx[ii],
				static_cast<float>(inputScalars.Ny[ii]) * inputScalars.dy[ii] + inputScalars.by[ii],
				static_cast<float>(inputScalars.Nz[ii]) * inputScalars.dz[ii] + inputScalars.bz[ii]);
		}
		if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
			erotusSens[0] = inputScalars.det_per_ring % local_size[0];
			erotusSens[1] = inputScalars.det_per_ring % local_size[1];
			if (erotusSens[1] > 0)
				erotusSens[1] = (local_size[1] - erotusSens[1]);
			if (erotusSens[0] > 0)
				erotusSens[0] = (local_size[0] - erotusSens[0]);
			d_xFull.resize(1);
			d_zFull.resize(1);
		}
		if (d_Summ.size() < 1)
			d_Summ.resize(1);
#if defined(OPENCL) || defined(METAL)
		region = { inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0] * inputScalars.nRekos };
#endif // END CUDA
		return 0;
		}

	/// <summary>
	/// This function first creates the necessary OpenCL/CUDA/HIP/Metal buffers and then writes the data to them
	/// </summary>
	/// <param name="length the number of measurements/projection/sinograms per subset"></param>
	/// <param name="x the x/y/z coordinates for the detectors (PET and SPECT) or source and detector (CT). z-coordinate applies only for CT"></param>
	/// <param name="z_det the z coordinates for the detectors (PET and SPECT) or the directional vectors for the detector panel pixels (CT)"></param>
	/// <param name="xy_index subset indices for subsets types &lt; 8, x/y dimensions"></param>
	/// <param name="z_index same as above but for z dimension"></param>
	/// <param name="L raw data detector indices"></param>
	/// <param name="pituus cumulative sum of length"></param>
	/// <param name="atten attenuation image"></param>
	/// <param name="norm normalization matrix"></param>
	/// <param name="extraCorr scatter data (for multiplicative scatter correction)"></param>
	/// <param name="x_center x-coordinates of the voxel centers"></param>
	/// <param name="y_center y-coordinates of the voxel centers"></param>
	/// <param name="z_center z-coordinates of the voxel centers"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="w_vec specifies some of the special options/parameters used"></param>
	/// <param name="MethodList specifies the algorithms and priors used"></param>
	/// <param name="Sin measurement data (sinograms or projections)"></param>
	/// <param name="sc_ra randoms and/or scatter data (for additive scatter correction or for randoms correction)"></param>
	/// <returns></returns>
	inline STATUS_t createAndWriteBuffers(const std::vector<int64_t>&length, const float* x, const float* z_det, const uint32_t * xy_index,
		const uint16_t * z_index, const uint16_t * L, const int64_t * pituus, const float* atten, const float* norm, const float* extraCorr,
		const scalarStruct & inputScalars, const Weighting & w_vec, const RecMethods & MethodList) {
		STATUS_t status = SUCCESS_VALUE;
		size_t vecSize = 1;
		if ((inputScalars.PET || inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0)
			vecSize = static_cast<size_t>(inputScalars.nRowsD) * static_cast<size_t>(inputScalars.nColsD);
		const size_t maskPriorDepth = inputScalars.multiResolution
			? static_cast<size_t>(inputScalars.Nz[0])
			: static_cast<size_t>(inputScalars.maskBPZ);
		// Check for the correct size of the inputs
		{
			bool missing = false;
			auto needNonEmpty = [&](const char* name, const size_t have) {
				if (have == 0) {
					mexPrintBase("%s is empty, but this configuration requires it\n", name);
					mexEval();
					missing = true;
				}
			};
			if (inputScalars.projector_type != 6) {
				if (inputScalars.BPType == 2 || inputScalars.BPType == 3 || inputScalars.FPType == 2 || inputScalars.FPType == 3)
					needNonEmpty("V (tube-of-response volume)", inputScalars.size_V);
				if ((!(inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0) || inputScalars.indexBased)
					needNonEmpty("x (detector coordinates)", inputScalars.size_of_x);
				if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
					needNonEmpty("x (full detector coordinates, required for the sensitivity image)", inputScalars.size_of_x);
					needNonEmpty("z (full detector coordinates, required for the sensitivity image)", inputScalars.size_z);
				}
			}
			for (uint32_t timestep = 0; timestep < inputScalars.Nt; timestep++)
				for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
					const size_t idx = static_cast<size_t>(kk) + static_cast<size_t>(timestep) * static_cast<size_t>(inputScalars.subsets);
					if (kk >= length.size() || idx >= length.size()) {
						mexPrintBase("length has %llu entries, but subset/timestep indexing needs at least ",
							static_cast<unsigned long long>(length.size()));
						mexPrintBase("%llu\n", static_cast<unsigned long long>(idx + 1));
						mexEval();
						missing = true;
					}
					else if (length[kk] <= 0 || length[idx] <= 0) {
						mexPrintBase("subset %u ", kk);
						mexPrintBase("of timestep %u contains no measurements\n", timestep);
						mexEval();
						missing = true;
					}
				}
			if (missing)
				return (STATUS_t)(-1);
		}
		// NLM anatomical reference image
		if (w_vec.NLM_anatomical && (MethodList.NLM || MethodList.ProxNLM)) {
			if (inputScalars.useImages) {
				CREATE_FLOAT_TEXTURE3D_FROM_HOST(d_urefIm, uRefArray, w_vec.NLM_ref, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0],
					BACKEND_TEXTURE_POINT, BACKEND_TEXTURE_DEFAULT_FLAGS);
			}
			else
				ALLOC_BUFFER(d_uref, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.im_dim[0]);
			CHECK(status, "\n", (STATUS_t)(-1));
			memAlloc.NLMRef = inputScalars.useImages ? 1 : 2;
		}
		// Create the prior image copy if needed, reuse the forward projection one otherwise
		// We define the image here and later input the current estimate to this image, if necessary
		if (inputScalars.projector_type == 6 || (inputScalars.FPType == 5 &&
			(MethodList.NLM || MethodList.RDP || MethodList.TV || MethodList.GGMRF || MethodList.APLS || MethodList.hyperbolic))) {
			if (inputScalars.useImages && !inputScalars.largeDim) {
				CREATE_FLOAT_TEXTURE3D_EMPTY(d_inputI, imArray, region[0], region[1], region[2]);
				CHECK(status, "Failed to create prior image\n", (STATUS_t)(-1));
			}
		}
		// RDP reference image
		if (MethodList.RDP && w_vec.RDPLargeNeighbor && w_vec.RDP_anatomical) {
			if (inputScalars.useImages) {
				CREATE_FLOAT_TEXTURE3D_EMPTY(d_RDPrefI, imArray, region[0], region[1], region[2]);
				CHECK(status, "Failed to create RDP reference image\n", (STATUS_t)(-1));
			}
		}
		// Create the necessary buffers
		// Distance-based weighting for GGMRF, RDP and hyperbolic prior
		if (MethodList.GGMRF || (MethodList.RDP && w_vec.RDPLargeNeighbor) || MethodList.hyperbolic) {
			ALLOC_BUFFER(d_weights, CL_MEM_READ_ONLY, sizeof(float) * (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1) - 1);
			CHECK(status, "\n", (STATUS_t)(-1));
			memAlloc.GGMRF = true;
		}
		memAlloc.tSteps = inputScalars.Nt;
		if ((inputScalars.useExtendedFOV && !inputScalars.multiResolution) || inputScalars.maskBP) {
			if (inputScalars.useBuffers) {
				ALLOC_BUFFER(d_maskPriorB, CL_MEM_READ_ONLY, sizeof(uint8_t) * static_cast<size_t>(inputScalars.Nx[0]) * static_cast<size_t>(inputScalars.Ny[0]) * maskPriorDepth);
			}
			else {
				if (inputScalars.maskBPZ > 1) {
					CREATE_MASK_TEXTURE3D_FROM_HOST(d_maskPrior, d_maskPrior3, maskArrayPrior, w_vec.maskPrior,
						inputScalars.Nx[0], inputScalars.Ny[0], maskPriorDepth, maskPriorDepth, maskPriorDepth, BACKEND_TEXTURE_READ_AS_INTEGER);
				}
				else {
					CREATE_MASK_TEXTURE2D_FROM_HOST(d_maskPrior, maskArrayPrior, w_vec.maskPrior,
						inputScalars.Nx[0], inputScalars.Ny[0], BACKEND_TEXTURE_READ_AS_INTEGER);
				}
				if (DEBUG) {
					mexPrintBase("imX = %u\n", inputScalars.Nx[0]);
					mexPrintBase("imY = %u\n", inputScalars.Ny[0]);
					mexPrintBase("imZ = %u\n", static_cast<unsigned int>(maskPriorDepth));
					mexEval();
				}
			}
			memAlloc.priorMask = true;
			CHECK(status, "\n", (STATUS_t)(-1));
		}
		if (inputScalars.projector_type != 6) {
			if (inputScalars.BPType == 2 || inputScalars.BPType == 3 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
				ALLOC_BUFFER(d_V, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.size_V);
				CHECK(status, "\n", (STATUS_t)(-1));
				memAlloc.V = true;
			}
			// Detector coordinates
			if ((!(inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0) || inputScalars.indexBased) {
				ALLOC_BUFFER(d_x[0][0], CL_MEM_READ_ONLY, sizeof(float) * inputScalars.size_of_x);
				CHECK(status, "\n", (STATUS_t)(-1));
				memAlloc.xSteps++;
			}
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				ALLOC_BUFFER(d_xFull[0], CL_MEM_READ_ONLY, sizeof(float) * inputScalars.size_of_x);
				CHECK(status, "\n", (STATUS_t)(-1));
				memAlloc.xFull = true;
			}
			// Mask images. Dynamic SPECT keeps separate resources per timestep;
			// detector-head stacks are shared within each timestep.
			if (inputScalars.maskFP || inputScalars.maskBP) {
				if (inputScalars.maskFP) {
					if (inputScalars.useBuffers) {
						d_maskFPB.resize(inputScalars.Nt);
						for (uint32_t timestep = 0; timestep < inputScalars.Nt; timestep++) {
							if (inputScalars.maskFPZ > 1 && !(inputScalars.SPECT && inputScalars.maskFPZ == inputScalars.nHeads)) {
								d_maskFPB[timestep].resize(inputScalars.subsetsUsed);
								for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
									const uint32_t indD = kk + timestep * inputScalars.subsets;
									ALLOC_BUFFER(d_maskFPB[timestep][kk], CL_MEM_READ_ONLY, sizeof(uint8_t) * inputScalars.nRowsD * inputScalars.nColsD * length[indD]);
									CHECK(status, "\n", (STATUS_t)(-1));
								}
							}
							else {
								d_maskFPB[timestep].resize(1);
								const size_t maskDepth = (inputScalars.SPECT && inputScalars.maskFPZ == inputScalars.nHeads) ? inputScalars.nHeads : 1ULL;
								ALLOC_BUFFER(d_maskFPB[timestep][0], CL_MEM_READ_ONLY, sizeof(uint8_t) * inputScalars.nRowsD * inputScalars.nColsD * maskDepth);
								CHECK(status, "\n", (STATUS_t)(-1));
							}
						}
					}
					else if (inputScalars.maskFPZ > 1) {
						d_maskFP3.resize(inputScalars.Nt);
						maskArrayFP.resize(inputScalars.Nt);
						for (uint32_t timestep = 0; timestep < inputScalars.Nt; timestep++) {
							if (inputScalars.SPECT && inputScalars.maskFPZ == inputScalars.nHeads) {
								d_maskFP3[timestep].resize(1);
								maskArrayFP[timestep].resize(1);
								CREATE_MASK_TEXTURE3D_FROM_HOST(d_maskFP3[timestep][0], d_maskFP3[timestep][0], maskArrayFP[timestep][0], w_vec.maskFP,
									inputScalars.nRowsD, inputScalars.nColsD, inputScalars.nHeads, inputScalars.nHeads, inputScalars.nHeads, BACKEND_TEXTURE_READ_AS_INTEGER);
								CHECK(status, "\n", (STATUS_t)(-1));
							}
							else {
								d_maskFP3[timestep].resize(inputScalars.subsetsUsed);
								maskArrayFP[timestep].resize(inputScalars.subsetsUsed);
								for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
									const uint32_t indD = kk + timestep * inputScalars.subsets;
									CREATE_MASK_TEXTURE3D_FROM_HOST(d_maskFP3[timestep][kk], d_maskFP3[timestep][kk], maskArrayFP[timestep][kk], &w_vec.maskFP[pituus[indD] * vecSize],
										inputScalars.nRowsD, inputScalars.nColsD, length[indD], length[indD], length[indD], BACKEND_TEXTURE_READ_AS_INTEGER);
									CHECK(status, "\n", (STATUS_t)(-1));
								}
							}
						}
					}
					else {
						maskArrayFP.resize(1);
						maskArrayFP[0].resize(1);
						CREATE_MASK_TEXTURE2D_FROM_HOST(d_maskFP, maskArrayFP[0][0], w_vec.maskFP,
							inputScalars.nRowsD, inputScalars.nColsD, BACKEND_TEXTURE_READ_AS_INTEGER);
						CHECK(status, "\n", (STATUS_t)(-1));
					}
					memAlloc.maskFP = true;
				}
				if (inputScalars.maskBP)
					memAlloc.maskBP = true;
			}
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				ALLOC_BUFFER(d_zFull[0], CL_MEM_READ_ONLY, sizeof(float) * inputScalars.size_z);
				CHECK(status, "\n", (STATUS_t)(-1));
				memAlloc.zFull = true;
			}
			if (inputScalars.SPECT) {
				ALLOC_BUFFER(d_rayShiftsDetector, CL_MEM_READ_ONLY, sizeof(float) * 2 * inputScalars.n_rays * inputScalars.n_rays3D * inputScalars.nRowsD *
					inputScalars.nColsD * inputScalars.nHeads);
				CHECK(status, "\n", (STATUS_t)(-1));
				ALLOC_BUFFER(d_rayShiftsSource, CL_MEM_READ_ONLY, sizeof(float) * 2 * inputScalars.n_rays * inputScalars.n_rays3D * inputScalars.nRowsD *
					inputScalars.nColsD * inputScalars.nHeads);
				CHECK(status, "\n", (STATUS_t)(-1));
				memAlloc.rayShifts = true;
			}
			if (inputScalars.eFOV && !inputScalars.multiResolution) {
				ALLOC_BUFFER(d_eFOVIndices, CL_MEM_READ_ONLY, sizeof(uint8_t) * inputScalars.Nz[0]);
				CHECK(status, "\n", (STATUS_t)(-1));
				memAlloc.eFOV = true;
			}
			if (inputScalars.CT && MethodList.FDK && inputScalars.useFDKWeights) {
				ALLOC_BUFFER(d_angle, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.nProjections);
				CHECK(status, "\n", (STATUS_t)(-1));
				memAlloc.angle = true;
			}
			// TOF bin centers
			if (inputScalars.TOF) {
				ALLOC_BUFFER(d_TOFCenter, CL_MEM_READ_ONLY, sizeof(float) * inputScalars.nBins);
				CHECK(status, "\n", (STATUS_t)(-1));
				memAlloc.TOF = true;
			}
			const auto maskBPFlags = (inputScalars.BPType == 4 && !inputScalars.CT) ? BACKEND_TEXTURE_NORMALIZED : BACKEND_TEXTURE_READ_AS_INTEGER;
			size_t maskBPTextureOffset = 0ULL;
			for (uint32_t timestep = 0; timestep < inputScalars.Nt; timestep++) {
				if (inputScalars.maskBP) {
					for (size_t volume = 0; volume < inputScalars.nMultiVolumes + 1; volume++) {
						const size_t maskBPDepth = inputScalars.maskBPZ > 1U
							? (inputScalars.multiResolution ? static_cast<size_t>(inputScalars.Nz[volume]) : static_cast<size_t>(inputScalars.maskBPZ))
							: 1ULL;
						const size_t maskBPElements = static_cast<size_t>(inputScalars.Nx[volume]) * static_cast<size_t>(inputScalars.Ny[volume]) * maskBPDepth;
						if (inputScalars.useBuffers) {
							ALLOC_BUFFER(d_maskBPB[timestep][volume], CL_MEM_READ_ONLY, sizeof(uint8_t) * maskBPElements);
						}
						else {
							if (inputScalars.maskBPZ > 1) {
								CREATE_MASK_TEXTURE3D_FROM_HOST(d_maskBP3[timestep][volume], d_maskBP3[timestep][volume], maskArrayBP[timestep][volume], &w_vec.maskBP[maskBPTextureOffset],
									inputScalars.Nx[volume], inputScalars.Ny[volume], maskBPDepth, maskBPDepth, maskBPDepth, maskBPFlags);
							}
							else {
								CREATE_MASK_TEXTURE2D_FROM_HOST(d_maskBP[timestep][volume], maskArrayBP[timestep][volume], &w_vec.maskBP[maskBPTextureOffset],
									inputScalars.Nx[volume], inputScalars.Ny[volume], maskBPFlags);
							}
							maskBPTextureOffset += maskBPElements;
						}
						CHECK(status, "\n", (STATUS_t)(-1));
					}
				}
				if (inputScalars.attenuation_correction && inputScalars.CTAttenuation) {
					if (inputScalars.size_atten > inputScalars.im_dim[0] || timestep == 0) {
						if (inputScalars.useBuffers)
							ALLOC_BUFFER(d_attenB[timestep], CL_MEM_READ_ONLY, sizeof(float) * inputScalars.im_dim[0]);
						else {
							const bool interpolationTexture = inputScalars.FPType == 4 || inputScalars.BPType == 4;
							CREATE_FLOAT_TEXTURE3D_FROM_HOST(d_attenIm[timestep], atArray, &atten[inputScalars.im_dim[0] * timestep],
								inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0],
								interpolationTexture ? BACKEND_TEXTURE_LINEAR : BACKEND_TEXTURE_POINT,
								interpolationTexture ? BACKEND_TEXTURE_NORMALIZED : BACKEND_TEXTURE_DEFAULT_FLAGS);
							CHECK(status, "\n", (STATUS_t)(-1));
						}
						memAlloc.atten = true;
						memAlloc.attenSize++;
					}
				}
				for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
					const uint32_t indD = kk + timestep * inputScalars.subsets;
					if (inputScalars.SPECT) {
						ALLOC_BUFFER(d_detectorVector[timestep][kk], CL_MEM_READ_ONLY, sizeof(uint32_t) * length[indD]);
						CHECK(status, "\n", (STATUS_t)(-1));
					}
					if (inputScalars.CT || inputScalars.SPECT) {
						ALLOC_BUFFER(d_x[timestep][kk], CL_MEM_READ_ONLY, sizeof(float) * length[indD] * 6);
						CHECK(status, "\n", (STATUS_t)(-1));
						memAlloc.xSteps++;
					}
					// First condition: load all data at once 
					// Second condition: load one subset at a time (only 1 buffer required, loadCoord reloads it for each subset/timestep)
					else if (inputScalars.listmode > 0 && !inputScalars.indexBased && (inputScalars.loadTOF || (kk == inputScalars.osa_iter0 && timestep == inputScalars.timestep0))) {
						ALLOC_BUFFER(d_x[timestep][kk], CL_MEM_READ_ONLY, sizeof(float) * length[kk + timestep * inputScalars.subsets] * 6);
						CHECK(status, "\n", (STATUS_t)(-1));
						memAlloc.xSteps++;
						if (DEBUG) {
							mexPrintBase("length[kk + timestep * inputScalars.subsets] * 6 = %u\n", length[kk + timestep * inputScalars.subsets] * 6);
							mexEval();
						}
					}
					if ((inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode != 1) {
						size_t coef = 2;
						if (inputScalars.useHelical)
							coef = 1;
						else if (inputScalars.pitch)
							coef = 6;
						ALLOC_BUFFER(d_z[timestep][kk], CL_MEM_READ_ONLY, sizeof(float) * length[indD] * coef);
						CHECK(status, "\n", (STATUS_t)(-1));
						memAlloc.zType = 1;
						memAlloc.zSteps++;
					}
					else {
						const uint32_t zTimestep = inputScalars.listmode > 0 ? 0U : timestep;
						const uint32_t zSubset = inputScalars.listmode > 0 ? 0U : kk;
						if (inputScalars.PET && inputScalars.listmode == 0) {
							if (inputScalars.nLayers > 1)
								ALLOC_BUFFER(d_z[timestep][kk], CL_MEM_READ_ONLY, sizeof(float) * length[kk] * 3);
							else
								ALLOC_BUFFER(d_z[timestep][kk], CL_MEM_READ_ONLY, sizeof(float) * length[kk] * 2);
							memAlloc.zType = 1;
							memAlloc.zSteps++;
						}
						else if (kk == inputScalars.osa_iter0 && (inputScalars.listmode == 0 || inputScalars.indexBased)) {
							ALLOC_BUFFER(d_z[zTimestep][zSubset], CL_MEM_READ_ONLY, sizeof(float) * (inputScalars.size_z > 0 ? inputScalars.size_z : static_cast<size_t>(1)));
							memAlloc.zType = 0;
							memAlloc.zSteps = kk;
						}
						else {
							ALLOC_BUFFER(d_z[zTimestep][zSubset], CL_MEM_READ_ONLY, sizeof(float) * (inputScalars.size_z > 0 ? inputScalars.size_z : static_cast<size_t>(1)));
							memAlloc.zType = 0;
							memAlloc.zSteps = kk;
						}
						CHECK(status, "\n", (STATUS_t)(-1));
					}
					if (inputScalars.size_scat > 1 && inputScalars.scatter == 1U) { // Scatter correction buffer
						ALLOC_BUFFER(d_scat[timestep][kk], CL_MEM_READ_ONLY, sizeof(float) * length[indD] * vecSize);
						CHECK(status, "\n", (STATUS_t)(-1));
						memAlloc.extra = true;
						memAlloc.eSteps++;
					}
					if (inputScalars.size_norm > 1 && inputScalars.normalization_correction) {
						const size_t normElements = (inputScalars.SPECT && inputScalars.normZ == inputScalars.nHeads)
							? static_cast<size_t>(inputScalars.nRowsD) * static_cast<size_t>(inputScalars.nColsD) * static_cast<size_t>(inputScalars.nHeads)
							: static_cast<size_t>(length[indD]) * vecSize;
						ALLOC_BUFFER(d_norm[timestep][kk], CL_MEM_READ_ONLY, sizeof(float) * normElements);
						CHECK(status, "\n", (STATUS_t)(-1));
						memAlloc.norm = true;
						memAlloc.nSteps++;
					}
					if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation) {
						ALLOC_BUFFER(d_atten[timestep][kk], CL_MEM_READ_ONLY, sizeof(float) * length[indD] * vecSize);
						CHECK(status, "\n", (STATUS_t)(-1));
						memAlloc.attenM = true;
						memAlloc.aSteps++;
					}
					if (inputScalars.listmode > 0 && inputScalars.indexBased) {
						if (inputScalars.loadTOF || (kk == inputScalars.osa_iter0 && !inputScalars.loadTOF && timestep == inputScalars.timestep0)) {
							ALLOC_BUFFER(d_trIndex[timestep][kk], CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk + timestep * inputScalars.subsets] * 2);
							CHECK(status, "\n", (STATUS_t)(-1));
							ALLOC_BUFFER(d_axIndex[timestep][kk], CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk + timestep * inputScalars.subsets] * 2);
							CHECK(status, "\n", (STATUS_t)(-1));
							memAlloc.indexBased = true;
							memAlloc.iSteps++;
						}
					}
					if (inputScalars.listmode > 0 && inputScalars.TOF) {
						if (inputScalars.loadTOF || (kk == inputScalars.osa_iter0 && !inputScalars.loadTOF && timestep == inputScalars.timestep0)) {
							ALLOC_BUFFER(d_TOFIndex[timestep][kk], CL_MEM_READ_ONLY, sizeof(uint8_t) * length[kk + timestep * inputScalars.subsets]);
							CHECK(status, "\n", (STATUS_t)(-1));
							memAlloc.TOFIndex = true;
							memAlloc.TOFSteps++;
						}
					}
				}
			}
			for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
				// Redundancy weighting
				if (inputScalars.offset && ((inputScalars.BPType == 4 && inputScalars.CT) || inputScalars.BPType == 5)) {
					ALLOC_BUFFER(d_T[kk], CL_MEM_READ_ONLY, sizeof(float) * length[kk]);
					CHECK(status, "\n", (STATUS_t)(-1));
					memAlloc.offsetT = true;
					memAlloc.oSteps++;
				}
#if !defined(METAL)
				if (inputScalars.BPType == 5 && inputScalars.CT && inputScalars.listmode == 0) {
					ALLOC_BUFFER(d_geomProj5[kk], CL_MEM_READ_ONLY, sizeof(float) * length[kk] * 16);
					CHECK(status, "\n", (STATUS_t)-1);
#if defined(CUDA) || defined(HIP)
					memAlloc.geom5 = true;
					memAlloc.g5Steps++;
#endif // END CUDA
				}
#endif
				// Indices corresponding to the detector index (Sinogram data) or the detector number (raw data) at each measurement
				// Note that raw data format is not used at the moment
				if (inputScalars.raw && inputScalars.listmode != 1) {
					ALLOC_BUFFER(d_L[kk], CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk] * 2);
					CHECK(status, "\n", (STATUS_t)(-1));
					memAlloc.raw = true;
					memAlloc.lSteps++;
				}
				else if (inputScalars.listmode != 1 && ((!inputScalars.CT && !inputScalars.SPECT && !inputScalars.PET) && 
					(inputScalars.subsets > 1 && (inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7)))) {
					ALLOC_BUFFER(d_xyindex[kk], CL_MEM_READ_ONLY, sizeof(uint32_t) * length[kk]);
					CHECK(status, "\n", (STATUS_t)(-1));
					ALLOC_BUFFER(d_zindex[kk], CL_MEM_READ_ONLY, sizeof(uint16_t) * length[kk]);
					CHECK(status, "\n", (STATUS_t)(-1));
					memAlloc.subInd = true;
					memAlloc.iSteps++;
				}
			}
		}
		CHECK(status, "Buffer creation failed\n", (STATUS_t)(-1));
		if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Buffer creation succeeded\n");
		}

		// assign values to the buffers
		if (MethodList.GGMRF || (MethodList.RDP && w_vec.RDPLargeNeighbor) || MethodList.hyperbolic) {
			WRITE_BUFFER(d_weights, sizeof(float) * (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1) - 1, w_vec.weights);
			memSize += (sizeof(float) * (w_vec.Ndx * 2 + 1) * (w_vec.Ndy * 2 + 1) * (w_vec.Ndz * 2 + 1));
			CHECK(status, "\n", (STATUS_t)(-1));
		}
		if (w_vec.NLM_anatomical && (MethodList.NLM || MethodList.ProxNLM)) {
			if (!inputScalars.useImages) {
				WRITE_BUFFER(d_uref, sizeof(float) * inputScalars.im_dim[0], w_vec.NLM_ref);
				CHECK(status, "\n", (STATUS_t)(-1));
			}
			memSize += (sizeof(float) * inputScalars.im_dim[0]);
		}
		if (inputScalars.projector_type != 6) {
			if (inputScalars.BPType == 2 || inputScalars.BPType == 3 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
				WRITE_BUFFER(d_V, sizeof(float) * inputScalars.size_V, inputScalars.V);
				CHECK(status, "\n", (STATUS_t)(-1));
				memSize += (sizeof(float) * inputScalars.size_V);
			}
			if ((!(inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0) || inputScalars.indexBased) {
				WRITE_BUFFER(d_x[0][0], sizeof(float) * inputScalars.size_of_x, x);
				CHECK(status, "\n", (STATUS_t)(-1));
				memSize += (sizeof(float) * inputScalars.size_of_x);
			}
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				WRITE_BUFFER(d_xFull[0], sizeof(float) * inputScalars.size_of_x, x);
				CHECK(status, "\n", (STATUS_t)(-1));
				WRITE_BUFFER(d_zFull[0], sizeof(float) * inputScalars.size_z, z_det);
				CHECK(status, "\n", (STATUS_t)(-1));
				memSize += (sizeof(float) * inputScalars.size_of_x + sizeof(float) * inputScalars.size_z);
			}
			if (inputScalars.eFOV && !inputScalars.multiResolution) {
				WRITE_BUFFER(d_eFOVIndices, sizeof(uint8_t) * inputScalars.Nz[0], w_vec.eFOVIndices);
				CHECK(status, "\n", (STATUS_t)(-1));
				memSize += (sizeof(uint8_t) * inputScalars.Nz[0]);
			}
			if (inputScalars.maskFP || inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution)) {
				if (inputScalars.useBuffers) {
					if (inputScalars.maskFP) {
						if (inputScalars.maskFPZ > 1 && !(inputScalars.SPECT && inputScalars.maskFPZ == inputScalars.nHeads))
							for (uint32_t timestep = 0; timestep < inputScalars.Nt; timestep++)
								for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
									const uint32_t indD = kk + timestep * inputScalars.subsets;
									WRITE_BUFFER(d_maskFPB[timestep][kk], sizeof(uint8_t) * inputScalars.nRowsD * inputScalars.nColsD * length[indD], &w_vec.maskFP[pituus[indD] * vecSize]);
								}
						else {
							for (uint32_t timestep = 0; timestep < inputScalars.Nt; timestep++)
								WRITE_BUFFER(d_maskFPB[timestep][0], sizeof(uint8_t) * inputScalars.nRowsD * inputScalars.nColsD *
									((inputScalars.SPECT && inputScalars.maskFPZ == inputScalars.nHeads) ? inputScalars.nHeads : 1ULL), w_vec.maskFP);
						}
						CHECK(status, "\n", (STATUS_t)(-1));
						const size_t maskFPDepth = inputScalars.maskFPZ > 1 ? static_cast<size_t>(inputScalars.maskFPZ) : 1ULL;
						memSize += sizeof(uint8_t) * static_cast<size_t>(inputScalars.nRowsD) * static_cast<size_t>(inputScalars.nColsD) * maskFPDepth;
					}
					if ((inputScalars.useExtendedFOV && !inputScalars.multiResolution) || inputScalars.maskBP) {
						WRITE_BUFFER(d_maskPriorB, sizeof(uint8_t) * static_cast<size_t>(inputScalars.Nx[0]) * static_cast<size_t>(inputScalars.Ny[0]) * maskPriorDepth, w_vec.maskPrior);
						CHECK(status, "\n", (STATUS_t)(-1));
						memSize += sizeof(uint8_t) * static_cast<size_t>(inputScalars.Nx[0]) * static_cast<size_t>(inputScalars.Ny[0]) * maskPriorDepth;
					}
				}
				else {
					if (inputScalars.maskFP) {
						const size_t maskFPDepth = inputScalars.maskFPZ > 1 ? static_cast<size_t>(inputScalars.maskFPZ) : 1ULL;
						memSize += (sizeof(uint8_t) * static_cast<size_t>(inputScalars.nRowsD) * static_cast<size_t>(inputScalars.nColsD) * maskFPDepth);
					}
					if ((inputScalars.useExtendedFOV && !inputScalars.multiResolution) || inputScalars.maskBP) {
						memSize += sizeof(uint8_t) * static_cast<size_t>(inputScalars.Nx[0]) * static_cast<size_t>(inputScalars.Ny[0]) * maskPriorDepth;
					}
				}
			}
			if (inputScalars.CT && MethodList.FDK && inputScalars.useFDKWeights) {
				WRITE_BUFFER(d_angle, sizeof(float) * inputScalars.nProjections, w_vec.angles);
				CHECK(status, "\n", (STATUS_t)(-1));
				memSize += (sizeof(float) * inputScalars.nProjections);
			}
			if (inputScalars.TOF) {
				WRITE_BUFFER(d_TOFCenter, sizeof(float) * inputScalars.nBins, inputScalars.TOFCenter);
				CHECK(status, "\n", (STATUS_t)(-1));
				memSize += (sizeof(float) * inputScalars.nBins);
			}
			if (inputScalars.SPECT) {
				const size_t rayShiftSize = static_cast<size_t>(2) * inputScalars.n_rays * inputScalars.n_rays3D * inputScalars.nRowsD * inputScalars.nColsD * inputScalars.nHeads;
				WRITE_BUFFER(d_rayShiftsDetector, sizeof(float) * rayShiftSize, w_vec.rayShiftsDetector);
				CHECK(status, "\n", (STATUS_t)(-1));
				memSize += sizeof(float) * rayShiftSize;
				WRITE_BUFFER(d_rayShiftsSource, sizeof(float) * rayShiftSize, w_vec.rayShiftsSource);
				CHECK(status, "\n", (STATUS_t)(-1));
				memSize += sizeof(float) * rayShiftSize;
			}

			if (DEBUG) {
				mexPrint("Timestep phase\n");
			}
			size_t maskBPBufferOffset = 0ULL;
			for (uint32_t timestep = 0; timestep < inputScalars.Nt; timestep++) {
				if (inputScalars.maskBP) {
					for (size_t volume = 0; volume < inputScalars.nMultiVolumes + 1; volume++) {
						const size_t maskBPDepth = inputScalars.maskBPZ > 1U
							? (inputScalars.multiResolution ? static_cast<size_t>(inputScalars.Nz[volume]) : static_cast<size_t>(inputScalars.maskBPZ))
							: 1ULL;
						const size_t maskBPElements = static_cast<size_t>(inputScalars.Nx[volume]) * static_cast<size_t>(inputScalars.Ny[volume]) * maskBPDepth;
						if (inputScalars.useBuffers) {
							WRITE_BUFFER(d_maskBPB[timestep][volume], sizeof(uint8_t) * maskBPElements, &w_vec.maskBP[maskBPBufferOffset]);
							CHECK(status, "\n", (STATUS_t)(-1));
							maskBPBufferOffset += maskBPElements;
						}
						memSize += (sizeof(uint8_t) * maskBPElements) / 1048576ULL;
					}
				}
				if (inputScalars.attenuation_correction && inputScalars.CTAttenuation) {
					if (inputScalars.useBuffers) {
						if (inputScalars.size_atten > inputScalars.im_dim[0]) {
							WRITE_BUFFER(d_attenB[timestep], sizeof(float) * inputScalars.im_dim[0], &atten[inputScalars.im_dim[0] * timestep]);
							CHECK(status, "\n", (STATUS_t)(-1));
							memSize += (sizeof(float) * inputScalars.im_dim[0]);
						}
						else if (timestep == 0) {
							WRITE_BUFFER(d_attenB[timestep], sizeof(float) * inputScalars.im_dim[0], atten);
							CHECK(status, "\n", (STATUS_t)(-1));
							memSize += (sizeof(float) * inputScalars.im_dim[0]);
						}
					}
				}
				for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
					const uint32_t indD = kk + timestep * inputScalars.subsets;
					if (inputScalars.SPECT) {
						WRITE_BUFFER(d_detectorVector[timestep][kk], sizeof(uint32_t) * length[indD],
							&w_vec.detectorVector[pituus[indD]]);
						CHECK(status, "\n", (STATUS_t)(-1));
						memSize += sizeof(uint32_t) * length[indD];
					}
					if ((inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0) {
						size_t kerroin = 2;
						if (inputScalars.pitch)
							kerroin = 6;
						else if (inputScalars.useHelical)
							kerroin = 1;
						WRITE_BUFFER(d_z[timestep][kk], sizeof(float) * length[indD] * kerroin, &z_det[pituus[indD] * kerroin]);
						CHECK(status, "\n", (STATUS_t)(-1));
						memSize += sizeof(float) * length[indD] * kerroin;
					}
					else {
						const uint32_t zTimestep = inputScalars.listmode > 0 ? 0U : timestep;
						const uint32_t zSubset = inputScalars.listmode > 0 ? 0U : kk;
						if (inputScalars.PET && inputScalars.listmode == 0) {
							int64_t kerroin = 2;
							if (inputScalars.nLayers > 1)
								kerroin = 3;
							WRITE_BUFFER(d_z[timestep][kk], sizeof(float) * length[kk] * kerroin, &z_det[pituus[kk] * kerroin]);
							memSize += (sizeof(float) * length[kk] * kerroin);
						}
						else if (kk == inputScalars.osa_iter0 && inputScalars.size_z > 0) {
							WRITE_BUFFER(d_z[zTimestep][zSubset], sizeof(float) * inputScalars.size_z, z_det);
							memSize += (sizeof(float) * inputScalars.size_z);
						}
						CHECK(status, "\n", (STATUS_t)(-1));
					}
					if ((inputScalars.CT || inputScalars.SPECT) && inputScalars.listmode == 0) {
						WRITE_BUFFER(d_x[timestep][kk], sizeof(float) * length[indD] * 6, &x[pituus[indD] * 6]);
						CHECK(status, "\n", (STATUS_t)(-1));
						memSize += sizeof(float) * length[indD] * 6;
					}
					else if (inputScalars.listmode > 0 && !inputScalars.indexBased) {
						if (inputScalars.loadTOF || (kk == inputScalars.osa_iter0 && timestep == inputScalars.timestep0)) {
							WRITE_BUFFER(d_x[timestep][kk], sizeof(float) * length[kk + timestep * inputScalars.subsets] * 6, 
								&w_vec.listCoord[pituus[kk + timestep * inputScalars.subsets] * 6]);
							CHECK(status, "\n", (STATUS_t)(-1));
							if (DEBUG) {
								mexPrintBase("length[kk + timestep * inputScalars.subsets] * 6 = %u\n", length[kk + timestep * inputScalars.subsets] * 6);
								mexPrintBase("pituus[kk + timestep * inputScalars.subsets] * 6 = %u\n", pituus[kk + timestep * inputScalars.subsets] * 6);
								mexPrintBase("w_vec.listCoord[pituus[kk + timestep * inputScalars.subsets] * 6] = %f\n", w_vec.listCoord[pituus[kk + timestep * inputScalars.subsets] * 6]);
								mexEval();
							}
							memSize += (sizeof(float) * length[kk] * 6);
						}
					}
					if (inputScalars.size_scat > 1ULL && inputScalars.scatter == 1U) { // Load scatter data
						WRITE_BUFFER(d_scat[timestep][kk], sizeof(float) * length[indD] * vecSize, &extraCorr[pituus[indD] * vecSize]);
						CHECK(status, "\n", (STATUS_t)(-1));
						memSize += sizeof(float) * length[indD] * vecSize;
					}
					if (inputScalars.size_norm > 1ULL && inputScalars.normalization_correction) {
						if (inputScalars.SPECT && inputScalars.normZ == inputScalars.nHeads)
							WRITE_BUFFER(d_norm[timestep][kk], sizeof(float) * inputScalars.nRowsD * inputScalars.nColsD * inputScalars.nHeads, norm);
						else
							WRITE_BUFFER(d_norm[timestep][kk], sizeof(float) * length[indD] * vecSize, &norm[pituus[indD] * vecSize]);
						CHECK(status, "\n", (STATUS_t)(-1));
						memSize += sizeof(float) * ((inputScalars.SPECT && inputScalars.normZ == inputScalars.nHeads)
							? static_cast<size_t>(inputScalars.nRowsD) * static_cast<size_t>(inputScalars.nColsD) * static_cast<size_t>(inputScalars.nHeads)
							: static_cast<size_t>(length[indD]) * vecSize);
					}
					if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation) {
						WRITE_BUFFER(d_atten[timestep][kk], sizeof(float) * length[indD] * vecSize, &atten[pituus[indD] * vecSize]);
						CHECK(status, "\n", (STATUS_t)(-1));
						memSize += sizeof(float) * length[indD] * vecSize;
					}
					if (inputScalars.listmode > 0 && inputScalars.indexBased) {
						// First condition: load all data at once. Second condition: load one subset at a time (only 1 buffer required for each timestep).
						if (inputScalars.loadTOF || (kk == inputScalars.osa_iter0 && !inputScalars.loadTOF && timestep == inputScalars.timestep0)) {
							WRITE_BUFFER(d_trIndex[timestep][kk], sizeof(uint16_t) * length[kk + timestep * inputScalars.subsets] * 2, 
								&w_vec.trIndex[pituus[kk + timestep * inputScalars.subsets] * 2]);
							CHECK(status, "\n", (STATUS_t)(-1));
							memSize += (sizeof(uint16_t) * length[kk] * 2);
							WRITE_BUFFER(d_axIndex[timestep][kk], sizeof(uint16_t) * length[kk + timestep * inputScalars.subsets] * 2, 
								&w_vec.axIndex[pituus[kk + timestep * inputScalars.subsets] * 2]);
							CHECK(status, "\n", (STATUS_t)(-1));
							memSize += (sizeof(uint16_t) * length[kk] * 2);
						}
					}
					if (inputScalars.listmode > 0 && inputScalars.TOF) {
						if (inputScalars.loadTOF || (kk == inputScalars.osa_iter0 && !inputScalars.loadTOF && timestep == inputScalars.timestep0)) {
							WRITE_BUFFER(d_TOFIndex[timestep][kk], sizeof(uint8_t) * length[kk + timestep * inputScalars.subsets], 
								&w_vec.TOFIndices[pituus[kk + timestep * inputScalars.subsets]]);
							CHECK(status, "\n", (STATUS_t)(-1));
							memSize += (sizeof(uint8_t) * length[kk]);
						}
					}
				}
			}
			for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
				if (inputScalars.offset && ((inputScalars.BPType == 4 && inputScalars.CT) || inputScalars.BPType == 5)) {
					WRITE_BUFFER(d_T[kk], sizeof(float) * length[kk], &inputScalars.T[pituus[kk]]);
					CHECK(status, "\n", (STATUS_t)(-1));
				}
				// Per projection: s (3), d3 (3), normX (3), normY (3), crossP (3), upperPart (1), i.e. 16 floats.
				// Previously computed in the kernel itself
#if !defined(METAL)
				if (inputScalars.BPType == 5 && inputScalars.CT && inputScalars.listmode == 0) {
					const float indX = static_cast<float>(inputScalars.nRowsD) / 2.f;
					const float indY = static_cast<float>(inputScalars.nColsD) / 2.f;
					const size_t uvStride = inputScalars.pitch ? 6 : 2;
					const float* xs = &x[pituus[kk] * 6];
					const float* uv = &z_det[pituus[kk] * uvStride];
					geomProj5Host[kk].resize(static_cast<size_t>(length[kk]) * 16);
					for (int64_t pp = 0; pp < length[kk]; pp++) {
						const float sX = xs[pp * 6], sY = xs[pp * 6 + 1], sZ = xs[pp * 6 + 2];
						const float dX = xs[pp * 6 + 3], dY = xs[pp * 6 + 4], dZ = xs[pp * 6 + 5];
						float aXx, aXy, aXz, aYx, aYy, aYz;
						if (inputScalars.pitch) {
							aXx = uv[pp * 6] * indX; aXy = uv[pp * 6 + 1] * indX; aXz = uv[pp * 6 + 2] * indX;
							aYx = uv[pp * 6 + 3] * indY; aYy = uv[pp * 6 + 4] * indY; aYz = uv[pp * 6 + 5] * indY;
						}
						else {
							aXx = uv[pp * 2] * indX; aXy = uv[pp * 2 + 1] * indX; aXz = 0.f;
							aYx = 0.f; aYy = 0.f; aYz = indY * w_vec.dPitchY;
						}
						// d3 = d - apuX - apuY
						const float d3x = dX - aXx - aYx, d3y = dY - aXy - aYy, d3z = dZ - aXz - aYz;
						// d2 = apuX - apuY
						const float d2x = aXx - aYx, d2y = aXy - aYy, d2z = aXz - aYz;
						const float nXl = std::sqrt(aXx * aXx + aXy * aXy + aXz * aXz);
						const float nYl = std::sqrt(aYx * aYx + aYy * aYy + aYz * aYz);
						// d3 - d = -apuX - apuY
						const float ddx = d3x - dX, ddy = d3y - dY, ddz = d3z - dZ;
						// crossP = cross(d2, d3 - d)
						const float cx = d2y * ddz - d2z * ddy;
						const float cy = d2z * ddx - d2x * ddz;
						const float cz = d2x * ddy - d2y * ddx;
						// upperPart = dot(crossP, s - d)
						const float up = cx * (sX - dX) + cy * (sY - dY) + cz * (sZ - dZ);
						float* g = &geomProj5Host[kk][pp * 16];
						g[0] = sX; g[1] = sY; g[2] = sZ;
						g[3] = d3x; g[4] = d3y; g[5] = d3z;
						g[6] = aXx / nXl; g[7] = aXy / nXl; g[8] = aXz / nXl;
						g[9] = aYx / nYl; g[10] = aYy / nYl; g[11] = aYz / nYl;
						g[12] = cx; g[13] = cy; g[14] = cz;
						g[15] = up;
					}
					WRITE_BUFFER(d_geomProj5[kk], sizeof(float) * length[kk] * 16, geomProj5Host[kk].data());
					CHECK(status, "\n", (STATUS_t)-1);
					memSize += (sizeof(float) * length[kk] * 16);
				}
#endif
				if (inputScalars.raw && inputScalars.listmode != 1) {
					WRITE_BUFFER(d_L[kk], sizeof(uint16_t) * length[kk] * 2, &L[pituus[kk] * 2]);
					CHECK(status, "\n", (STATUS_t)(-1));
					memSize += (sizeof(uint16_t) * length[kk] * 2);
				}
				else if (inputScalars.listmode != 1 && ((!inputScalars.CT && !inputScalars.SPECT && !inputScalars.PET) && (inputScalars.subsets > 1 && 
					(inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7)))) {
					WRITE_BUFFER(d_zindex[kk], sizeof(uint16_t) * length[kk], &z_index[pituus[kk]]);
					CHECK(status, "\n", (STATUS_t)(-1));
					WRITE_BUFFER(d_xyindex[kk], sizeof(uint32_t) * length[kk], &xy_index[pituus[kk]]);
					CHECK(status, "\n", (STATUS_t)(-1));
					memSize += (sizeof(uint32_t) * length[kk] + sizeof(uint16_t) * length[kk]);
				}
			}
			FINISH_QUEUE(status, "Buffer write failed\n", (STATUS_t)(-1));
		}
		CHECK(status, "Buffer write failed\n", (STATUS_t)(-1));
		if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Buffer write succeeded\n");
		}
		return SUCCESS_VALUE;
	}

	/// <summary>
	/// Resizes required vectors and then calls the function to create and write buffers. Also creates two necessary images
	/// </summary>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="w_vec specifies some of the special options/parameters used"></param>
	/// <param name="x the x/y/z coordinates for the detectors (PET and SPECT) or source and detector (CT). z-coordinate applies only for CT"></param>
	/// <param name="z_det the z coordinates for the detectors (PET and SPECT) or the directional vectors for the detector panel pixels (CT)"></param>
	/// <param name="xy_index subset indices for subsets types &lt; 8, x/y dimensions"></param>
	/// <param name="z_index same as above but for z dimension"></param>
	/// <param name="lor1 LORs to be discarded, i.e. the ray does not intersect the image"></param>
	/// <param name="L raw data detector indices"></param>
	/// <param name="pituus cumulative sum of length"></param>
	/// <param name="atten attenuation image"></param>
	/// <param name="norm normalization matrix"></param>
	/// <param name="extraCorr scatter data (for multiplicative scatter correction)"></param>
	/// <param name="V precomputed volume values for the volume of intersection based projector"></param>
	/// <param name="x_center x-coordinates of the voxel centers"></param>
	/// <param name="y_center y-coordinates of the voxel centers"></param>
	/// <param name="z_center z-coordinates of the voxel centers"></param>
	/// <param name="sc_ra randoms and/or scatter data (for additive scatter correction or for randoms correction)"></param>
	/// <param name="TOFCenter TOF bin center values"></param>
	/// <param name="length the number of measurements/projection/sinograms per subset"></param>
	/// <param name="Sin measurement data (sinograms or projections)"></param>
	/// <param name="reko_type for reconstruction algorithms requiring unique operations in FP or BP"></param>
	/// <param name="MethodList specifies the algorithms and priors used"></param>
	/// <returns></returns>
	inline int createBuffers(scalarStruct & inputScalars, Weighting & w_vec, const float* x, const float* z_det, const uint32_t * xy_index,
		const uint16_t * z_index, const uint16_t * L, const int64_t * pituus, const float* atten, const float* norm, const float* extraCorr,
		const std::vector<int64_t>&length, const RecMethods & MethodList, const int type = 0) {
		STATUS_t status = SUCCESS_VALUE;
#if defined(METAL)
		if (type == 0 && d_Summ.empty()) {
			d_Summ.resize(inputScalars.nMultiVolumes + 1);
			float zero = 0.f;
			for (uint32_t ii = 0; ii <= inputScalars.nMultiVolumes; ++ii) {
				ALLOC_BUFFER(d_Summ[ii], CL_MEM_READ_WRITE, sizeof(zero));
				CHECK(status, "Failed to create Metal sensitivity placeholder\n", -1);
				WRITE_BUFFER(d_Summ[ii], sizeof(zero), &zero);
				CHECK(status, "Failed to initialize Metal sensitivity placeholder\n", -1);
			}
		}
#elif defined(OPENCL)
		// The sensitivity image is only allocated (through ArrayFire) when it is actually computed
		// On GPUs, it doesn't matter if d_Summ is not allocated if it's not used, but with CPUs
		// it needs to be allocated
		// This affects all algorithms that don't specifically compute the sensitivity image
		if (d_Summ.size() < static_cast<size_t>(inputScalars.nMultiVolumes) + 1ULL)
			d_Summ.resize(inputScalars.nMultiVolumes + 1);
		for (uint32_t ii = 0; ii <= inputScalars.nMultiVolumes; ++ii) {
			// Only the entries that are not already backed by an ArrayFire allocation
			if (d_Summ[ii]() == NULL) {
				ALLOC_BUFFER(d_Summ[ii], CL_MEM_READ_WRITE, sizeof(float));
				OCL_CHECK(status, "Failed to create the sensitivity image placeholder\n", -1);
			}
		}
#endif // END METAL
		if (inputScalars.raw)
			d_L.resize(inputScalars.subsetsUsed);
		if ((inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7) && inputScalars.subsets > 1) {
			d_xyindex.resize(inputScalars.subsetsUsed);
			d_zindex.resize(inputScalars.subsetsUsed);
		}
			if (inputScalars.normalization_correction)
				d_norm.resize(inputScalars.Nt);
			if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation)
				d_atten.resize(inputScalars.Nt);
			if (inputScalars.attenuation_correction && inputScalars.CTAttenuation) {
				d_attenB.resize(inputScalars.Nt);
				d_attenIm.resize(inputScalars.Nt);
			}
			if (inputScalars.maskBP) {
				if (inputScalars.useBuffers)
					d_maskBPB.resize(inputScalars.Nt);
				else {
					maskArrayBP.resize(inputScalars.Nt);
					if (inputScalars.maskBPZ > 1)
						d_maskBP3.resize(inputScalars.Nt);
					else
						d_maskBP.resize(inputScalars.Nt);
				}
				for (int tt = 0; tt < inputScalars.Nt; tt++) {
					if (inputScalars.useBuffers)
						d_maskBPB[tt].resize(inputScalars.nMultiVolumes + 1);
					else {
						maskArrayBP[tt].resize(inputScalars.nMultiVolumes + 1);
						if (inputScalars.maskBPZ > 1)
							d_maskBP3[tt].resize(inputScalars.nMultiVolumes + 1);
						else
							d_maskBP[tt].resize(inputScalars.nMultiVolumes + 1);
					}
				}
			}
			for (int tt = 0; tt < inputScalars.Nt; tt++) {
				if (inputScalars.normalization_correction)
					d_norm[tt].resize(inputScalars.subsetsUsed);
				if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation)
					d_atten[tt].resize(inputScalars.subsetsUsed);
			}
			if (inputScalars.projector_type != 6) {
			d_scat.resize(inputScalars.Nt);
			d_x.resize(inputScalars.Nt);
			d_z.resize(inputScalars.Nt);
            d_detectorVector.resize(inputScalars.Nt);
			d_trIndex.resize(inputScalars.Nt);
			d_axIndex.resize(inputScalars.Nt);
				d_TOFIndex.resize(inputScalars.Nt);
				for (int tt = 0; tt < inputScalars.Nt; tt++) {
					d_scat[tt].resize(inputScalars.subsetsUsed);
				d_x[tt].resize(inputScalars.subsetsUsed);
				d_z[tt].resize(inputScalars.subsetsUsed);
                d_detectorVector[tt].resize(inputScalars.subsetsUsed);
				d_trIndex[tt].resize(inputScalars.subsetsUsed);
				d_axIndex[tt].resize(inputScalars.subsetsUsed);
				d_TOFIndex[tt].resize(inputScalars.subsetsUsed);
			}
		}
		if (inputScalars.offset && ((inputScalars.BPType == 4 && inputScalars.CT) || inputScalars.BPType == 5))
			d_T.resize(inputScalars.subsetsUsed);
#if !defined(METAL)
		if (inputScalars.BPType == 5 && inputScalars.CT && inputScalars.listmode == 0) {
			d_geomProj5.resize(inputScalars.subsetsUsed);
			geomProj5Host.resize(inputScalars.subsetsUsed);
		}
#endif

		status = createAndWriteBuffers(length, x, z_det, xy_index, z_index, L, pituus, atten, norm, extraCorr, inputScalars, w_vec, MethodList);
		if (status != SUCCESS_VALUE) {
			return status;
		}
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrintBase("Allocated a total of %u MB\n", memSize / 1048576ULL);
		return 0;
	}

	/// <summary>
	/// Inputs constant values to the kernels, i.e. values that do not change in each time step or iteration
	/// </summary>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="w_vec specifies some of the special options/parameters used"></param>
	/// <returns></returns>
	inline int initializeKernel(scalarStruct & inputScalars, Weighting & w_vec) {
#if defined(METAL)
		kParams.nRowsD = inputScalars.nRowsD;
		kParams.nColsD = inputScalars.nColsD;
		kParams.dPitch = dPitch;
		kParams.ellipseCenter = make_vec3<FLOAT3_t>(inputScalars.ellipseCenterX, inputScalars.ellipseCenterY, inputScalars.ellipseCenterZ);
		kParams.ellipseRadii = make_vec3<FLOAT3_t>(inputScalars.ellipseRadiusX, inputScalars.ellipseRadiusY, inputScalars.ellipseRadiusZ);
		kParams.ellipsePower = inputScalars.ellipsePower;
		kParams.dL = inputScalars.dL;
		kParams.global_factor = inputScalars.global_factor;
		kParams.epps = inputScalars.epps;
		kParams.det_per_ring = inputScalars.det_per_ring;
		kParams.sigma_x = inputScalars.sigma_x;
		kParams.coneOfResponseStdCoeffA = inputScalars.coneOfResponseStdCoeffA;
		kParams.coneOfResponseStdCoeffB = inputScalars.coneOfResponseStdCoeffB;
		kParams.coneOfResponseStdCoeffC = inputScalars.coneOfResponseStdCoeffC;
		kParams.tube_width = inputScalars.tube_width;
		kParams.cylRadiusProj3 = inputScalars.cylRadiusProj3;
		kParams.bmin = inputScalars.bmin;
		kParams.bmax = inputScalars.bmax;
		kParams.Vmax = inputScalars.Vmax;
		kParams.rings = inputScalars.rings;
		kParams.helicalRadius = inputScalars.helicalRadius;
		return 0;
#else
		// Set the kernelFP parameters that do not change
		if (inputScalars.FPType == 4 || inputScalars.FPType == 5 || inputScalars.FPType == 7) {
			KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.nRowsD);
			KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.nColsD);
			KARG(FPArgs, kernelFP, kernelIndFP, dPitch);
			if (inputScalars.useHelical) {
				KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.helicalRadius);
			}
		}

		if (inputScalars.BPType == 4 || inputScalars.BPType == 5 || inputScalars.BPType == 7) {
			KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.nRowsD);
			KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.nColsD);
			KARG(BPArgs, kernelBP, kernelIndBP, dPitch);
			if (inputScalars.useHelical)
				KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.helicalRadius);
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.nRowsD);
				KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.nColsD);
				KARG(SensArgs, kernelSensList, kernelIndSens, dPitch);
			}
		}
		if (inputScalars.FPType == 4) {
			KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.dL);
			KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.global_factor);
		}
		if (inputScalars.BPType == 4 && !inputScalars.CT) {
			KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.dL);
			KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.global_factor);
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.dL);
				KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.global_factor);
			}
		}
		if (inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
			KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.global_factor);
			KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.epps);
			KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.nRowsD);
			KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.det_per_ring);
			KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.sigma_x);
			if (inputScalars.SPECT) {
				KARG(FPArgs, kernelFP, kernelIndFP, d_rayShiftsDetector);
				KARG(FPArgs, kernelFP, kernelIndFP, d_rayShiftsSource);
				KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.coneOfResponseStdCoeffA);
				KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.coneOfResponseStdCoeffB);
				KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.coneOfResponseStdCoeffC);
				KARG(FPArgs, kernelFP, kernelIndFP, ellipseCenter);
				KARG(FPArgs, kernelFP, kernelIndFP, ellipseRadii);
				KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.ellipsePower);
			}
			KARG(FPArgs, kernelFP, kernelIndFP, dPitch);
			if (inputScalars.FPType == 2 || inputScalars.FPType == 3) {
				if (inputScalars.FPType == 2) {
					KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.tube_width);
				}
				else {
					KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.cylRadiusProj3);
				}
				KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.bmin);
				KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.bmax);
				KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.Vmax);
			}
		}
		if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
			KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.global_factor);
			KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.epps);
			KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.nRowsD);
			KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.det_per_ring);
			KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.sigma_x);
			if (inputScalars.SPECT) {
				KARG(BPArgs, kernelBP, kernelIndBP, d_rayShiftsDetector);
				KARG(BPArgs, kernelBP, kernelIndBP, d_rayShiftsSource);
				KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.coneOfResponseStdCoeffA);
				KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.coneOfResponseStdCoeffB);
				KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.coneOfResponseStdCoeffC);
				KARG(BPArgs, kernelBP, kernelIndBP, ellipseCenter);
				KARG(BPArgs, kernelBP, kernelIndBP, ellipseRadii);
				KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.ellipsePower);
			}
			KARG(BPArgs, kernelBP, kernelIndBP, dPitch);
			if (inputScalars.BPType == 2 || inputScalars.BPType == 3) {
				if (inputScalars.BPType == 2) {
					KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.tube_width);
				}
				else {
					KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.cylRadiusProj3);
				}
				KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.bmin);
				KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.bmax);
				KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.Vmax);
			}
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.global_factor);
				KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.epps);
				KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.nRowsD);
				KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.det_per_ring);
				KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.sigma_x);
				KARG(SensArgs, kernelSensList, kernelIndSens, dPitch);
				if (inputScalars.BPType == 2 || inputScalars.BPType == 3) {
					if (inputScalars.BPType == 2) {
						KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.tube_width);
					}
					else
						KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.cylRadiusProj3);
					KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.bmin);
					KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.bmax);
					KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.Vmax);
				}
			}
		}
		if (DEBUG) {
			mexPrintBase("inputScalars.nBins = %u\n", inputScalars.nBins);
			mexPrintBase("inputScalars.helicalRadius = %f\n", inputScalars.helicalRadius);
			mexEval();
		}
		if (inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3) {
			if (inputScalars.TOF) {
				KARG(FPArgs, kernelFP, kernelIndFP, d_TOFCenter);
			}
			if (inputScalars.FPType == 2 || inputScalars.FPType == 3) {
				KARG(FPArgs, kernelFP, kernelIndFP, d_V);
			}
			KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.nColsD);
		}
		if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
			if (inputScalars.TOF)
				KARG(BPArgs, kernelBP, kernelIndBP, d_TOFCenter);
			if (inputScalars.BPType == 2 || inputScalars.BPType == 3) {
				KARG(BPArgs, kernelBP, kernelIndBP, d_V);
			}
			KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.nColsD);
			if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
				if (inputScalars.TOF)
					KARG(SensArgs, kernelSensList, kernelIndSens, d_TOFCenter);
				if (inputScalars.BPType == 2 || inputScalars.BPType == 3) {
					KARG(SensArgs, kernelSensList, kernelIndSens, d_V);
				}
				KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.nColsD);
			}
		}
		if ((inputScalars.BPType == 4 || inputScalars.FPType == 4) && !inputScalars.CT && inputScalars.TOF) {
			if (inputScalars.FPType == 4) {
				KARG(FPArgs, kernelFP, kernelIndFP, d_TOFCenter);
				KARG(FPArgs, kernelFP, kernelIndFP, inputScalars.sigma_x);
			}
			if (inputScalars.BPType == 4) {
				KARG(BPArgs, kernelBP, kernelIndBP, d_TOFCenter);
				KARG(BPArgs, kernelBP, kernelIndBP, inputScalars.sigma_x);
				if (inputScalars.listmode > 0 && inputScalars.computeSensImag) {
					KARG(SensArgs, kernelSensList, kernelIndSens, d_TOFCenter);
					KARG(SensArgs, kernelSensList, kernelIndSens, inputScalars.sigma_x);
				}
			}
		}
		if (DEBUG) {
#if defined(CUDA) || defined(HIP)
			mexPrintBase("kernelIndFP = %u\n", FPArgs.size());
			mexPrintBase("kernelIndBP = %u\n", BPArgs.size());
#elif defined(OPENCL)
			mexPrintBase("kernelIndFP = %u\n", kernelIndFP);
			mexPrintBase("kernelIndBP = %u\n", kernelIndBP);
#endif // END CUDA
			mexEval();
		}
		return 0;
#endif
	}

	template <typename T>
	inline int loadCoord(uint32_t currentSubset, scalarStruct & inputScalars, const int64_t length, const T * listCoord, 
		const T * listCoordAx = nullptr, const uint8_t * TOFIndices = nullptr) {
		STATUS_t status = SUCCESS_VALUE;
		if (inputScalars.indexBased) {
#if defined(CUDA) || defined(HIP)
			// CUDA/HIP must explicitly release the previous allocation; OpenCL's cl::Buffer reassignment does this itself.
			getErrorString(cuMemFree(d_trIndex[0][0]));
			getErrorString(cuMemFree(d_axIndex[0][0]));
#endif // END CUDA
			ALLOC_BUFFER(d_trIndex[0][0], CL_MEM_READ_ONLY, sizeof(uint16_t) * length * 2);
			CHECK(status, "\n", -1);
			ALLOC_BUFFER(d_axIndex[0][0], CL_MEM_READ_ONLY, sizeof(uint16_t) * length * 2);
			CHECK(status, "\n", -1);
			WRITE_BUFFER(d_trIndex[0][0], sizeof(uint16_t) * length * 2, listCoord);
			CHECK(status, "\n", -1);
			WRITE_BUFFER(d_axIndex[0][0], sizeof(uint16_t) * length * 2, listCoordAx);
			CHECK(status, "\n", -1);
		}
		else {
#if defined(CUDA) || defined(HIP)
			getErrorString(cuMemFree(d_x[0][0]));
#endif // END CUDA
			ALLOC_BUFFER(d_x[0][0], CL_MEM_READ_ONLY, sizeof(float) * length * 6);
			CHECK(status, "\n", -1);
			WRITE_BUFFER(d_x[0][0], sizeof(float) * length * 6, listCoord);
			CHECK(status, "\n", -1);
		}
		if (inputScalars.TOF) {
#if defined(CUDA) || defined(HIP)
			getErrorString(cuMemFree(d_TOFIndex[0][0]));
#endif // END CUDA
			ALLOC_BUFFER(d_TOFIndex[0][0], CL_MEM_READ_ONLY, sizeof(uint8_t) * length);
			CHECK(status, "\n", -1);
			WRITE_BUFFER(d_TOFIndex[0][0], sizeof(uint8_t) * length, TOFIndices);
			CHECK(status, "\n", -1);
		}
		return 0;
#if defined(OPENCL)
	}

	int computeConvolutionF(const scalarStruct & inputScalars, const int ii = 0) {
		STATUS_t status = SUCCESS_VALUE;
		cl::NDRange	globalC = { inputScalars.Nx[ii] + erotusBP[0][ii], inputScalars.Ny[ii] + erotusBP[1][ii], inputScalars.Nz[ii] };

		FINISH_QUEUE(status, "\n", -1);
		UINT32_t kernelInd = 0U;
		KARG(kArgs, kernelPSFf, kernelInd, d_imFinal[ii]);
		KARG(kArgs, kernelPSFf, kernelInd, d_imTemp[ii]);
		KARG(kArgs, kernelPSFf, kernelInd, d_g);
		KARG(kArgs, kernelPSFf, kernelInd, inputScalars.g_dim_x);
		KARG(kArgs, kernelPSFf, kernelInd, inputScalars.g_dim_y);
		KARG(kArgs, kernelPSFf, kernelInd, inputScalars.g_dim_z);
		status = CLCommandQueue[0].enqueueNDRangeKernel(kernelPSFf, cl::NDRange(), globalC, localPrior, NULL);
		OCL_CHECK(status, "\n", -1);
		if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Convolution kernel launched successfully\n");
		}
		FINISH_QUEUE(status, "\n", -1);
		if (DEBUG || inputScalars.verbose >= 3)
			mexPrint("Convolution computed");

		return status;
	}

	template <typename T>
	int computeConvolution(const scalarStruct & inputScalars, cl::Buffer & input, const int ii = 0, const T & cType = 0) {
		STATUS_t status = SUCCESS_VALUE;
		cl::NDRange	globalC = { inputScalars.Nx[ii] + erotusBP[0][ii], inputScalars.Ny[ii] + erotusBP[1][ii], inputScalars.Nz[ii] };

		cl::Buffer d_BPApu = cl::Buffer(CLContext, CL_MEM_READ_WRITE, sizeof(T) * inputScalars.im_dim[ii], NULL, &status);
		FINISH_QUEUE(status, "\n", -1);
		UINT32_t kernelInd = 0U;
		KARG(kArgs, kernelPSF, kernelInd, input);
		KARG(kArgs, kernelPSF, kernelInd, d_BPApu);
		KARG(kArgs, kernelPSF, kernelInd, d_g);
		KARG(kArgs, kernelPSF, kernelInd, inputScalars.g_dim_x);
		KARG(kArgs, kernelPSF, kernelInd, inputScalars.g_dim_y);
		KARG(kArgs, kernelPSF, kernelInd, inputScalars.g_dim_z);
		status = CLCommandQueue[0].enqueueNDRangeKernel(kernelPSF, cl::NDRange(), globalC, localPrior, NULL);
		OCL_CHECK(status, "\n", -1);
		if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Convolution kernel launched successfully\n");
		}
		FINISH_QUEUE(status, "\n", -1);
		if (DEBUG || inputScalars.verbose >= 3)
			mexPrint("Convolution computed");
		status = CLCommandQueue[0].enqueueCopyBuffer(d_BPApu, input, 0, 0, sizeof(T) * inputScalars.im_dim[ii]);
		OCL_CHECK(status, "\n", -1);
		FINISH_QUEUE(status, "\n", -1);

		return status;
	}

	int computeForward(const scalarStruct & inputScalars, const std::vector<int64_t>&length, const uint32_t osa_iter, const uint32_t timestep = 0) {
		STATUS_t status = SUCCESS_VALUE;
		cl::NDRange localF = { 64, 1, 1 };
		int indD = osa_iter + timestep * inputScalars.subsets;
		if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
			global = { inputScalars.nRowsD * inputScalars.nColsD * static_cast<size_t>(length[indD]) * inputScalars.nBins, 1, 1 };
		else
			if (inputScalars.listmode == 0)
				global = { static_cast<cl::size_type>(length[indD]) * inputScalars.nBins, 1, 1 };
			else
				global = { static_cast<cl::size_type>(length[indD]), 1, 1 };

		size_t erotusF = global[0] % localF[0];
		if (erotusF > 0)
			erotusF = localF[0] - erotusF;
		global = { static_cast<cl::size_type>(global[0] + erotusF), 1, 1 };
		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", global[0]);
			mexPrintBase("global[1] = %u\n", global[1]);
			mexPrintBase("global[2] = %u\n", global[2]);
			mexPrintBase("erotus[0] = %u\n", erotus[0]);
			mexPrintBase("erotus[1] = %u\n", erotus[1]);
			mexPrintBase("global.dimensions() = %u\n", global.dimensions());
			mexPrintBase("length[indD] = %u\n", length[indD]);
			mexPrintBase("listmode = %u\n", inputScalars.listmode);
			mexEval();
		}
		UINT32_t kernelInd = 0U;
		KARG(kArgs, kernelForward, kernelInd, d_output);
		KARG(kArgs, kernelForward, kernelInd, d_meas[osa_iter]);
		if (inputScalars.CT)
			KARG(kArgs, kernelForward, kernelInd, d_outputCT);
		if (inputScalars.randoms_correction)
			KARG(kArgs, kernelForward, kernelInd, d_rand[osa_iter]);
		status = CLCommandQueue[0].enqueueNDRangeKernel(kernelForward, cl::NDRange(), global, localF, NULL);
		OCL_CHECK(status, "\n", -1);
		if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Forward step kernel launched successfully\n");
		}
		FINISH_QUEUE(status, "\n", -1);
		if (inputScalars.verbose >= 3)
			mexPrint("Forward step computed");

		return status;
	}

	int computeEstimate(const scalarStruct & inputScalars, const int ii = 0, const int uu = 0, const int timestep = 0) {
		STATUS_t status = SUCCESS_VALUE;
		global = { inputScalars.Nx[ii] + erotusBP[0][ii], inputScalars.Ny[ii] + erotusBP[1][ii], inputScalars.Nz[ii] };

		UINT32_t kernelInd = 0U;
		KARG(kArgs, kernelEstimate, kernelInd, d_Summ[uu]);
		KARG(kArgs, kernelEstimate, kernelInd, vec_opencl.d_rhs_os[ii]);
		KARG(kArgs, kernelEstimate, kernelInd, d_imFinal[ii]);
		KARG(kArgs, kernelEstimate, kernelInd, inputScalars.epps);
		KARG(kArgs, kernelEstimate, kernelInd, d_N[ii]);
		KARG(kArgs, kernelEstimate, kernelInd, no_norm);
		if (inputScalars.CT)
			KARG(kArgs, kernelEstimate, kernelInd, inputScalars.flat);
		status = CLCommandQueue[0].enqueueNDRangeKernel(kernelEstimate, cl::NDRange(), global, localPrior, NULL);
		OCL_CHECK(status, "\n", -1);
		if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Estimate step kernel launched successfully\n");
		}
		FINISH_QUEUE(status, "\n", -1);
		if (inputScalars.verbose >= 3)
			mexPrint("Estimate step computed");

		return status;
#endif // END CUDA
	}

	/// <summary>
	/// Compute the forward projection for the selected projector type
	/// </summary>
	/// <param name="vec image estimates and backprojection"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="w_vec specifies some of the special options/parameters used"></param>
	/// <param name="outputFP the output forward projection array"></param>
	/// <param name="osa_iter current subset (sub-iteration)"></param>
	/// <param name="length the number of measurements/projection/sinograms per subset"></param>
	/// <param name="m_size for projector types 1-3, the total number of LORs"></param>
	/// <returns></returns>
#if defined(CUDA) || defined(HIP)
	inline int forwardProjection(scalarStruct & inputScalars, Weighting & w_vec, uint32_t osa_iter, uint32_t timestep, 
		const std::vector<int64_t>&length1, uint64_t m_size, int ii = 0, const int uu = 0) {
#elif defined(METAL) || defined(OPENCL)
	inline int forwardProjection(const scalarStruct & inputScalars, Weighting & w_vec, const uint32_t osa_iter, const uint32_t timestep, 
		const std::vector<int64_t>&length, uint64_t m_size, const int32_t ii = 0, const int uu = 0) {
#endif // END CUDA
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrintVar("Starting forward projection for projector type = ", inputScalars.FPType);
		STATUS_t status = SUCCESS_VALUE;
#if defined(CUDA) || defined(HIP)
		std::vector<int64_t> length = length1;
		std::vector<void*> kTemp = FPArgs;
#elif defined(OPENCL)
		kernelIndFPSubIter = kernelIndFP;
#elif defined(METAL)
		if (!queueFP || !kernelFP) {
			mexPrint("Unable to create Metal forward-projection encoder");
			return -1;
		}
		NS::SharedPtr<MTL::CommandBuffer> commandBuffer = NS::RetainPtr(queueFP->commandBuffer());
		NS::SharedPtr<MTL::ComputeCommandEncoder> encoder = NS::RetainPtr(commandBuffer->computeCommandEncoder());
		if (!commandBuffer || !encoder) {
			mexPrint("Unable to create Metal forward-projection encoder");
			return -1;
		}
		encoder->setComputePipelineState(kernelFP.get());
#endif // END CUDA
		// Per-launch copy of the work-group range: the 1D branch below has to flatten it, and local is a
		// member shared with the backprojection and with the other, multidimensional launches here.
		WORKRANGE_t localFP = local;
		if (inputScalars.FPType == 5) {
			SET_LAUNCH_RANGE3(global, inputScalars.nRowsD + erotus[0], (inputScalars.nColsD + NVOXELSFP - 1) / NVOXELSFP + erotus[1],
				length[osa_iter + timestep * inputScalars.subsets], localFP);
		}
		else if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0) {
			SET_LAUNCH_RANGE3(global, inputScalars.nRowsD + erotus[0], inputScalars.nColsD + erotus[1],
				length[osa_iter + timestep * inputScalars.subsets], localFP);
		}
		else {
			erotus[0] = length[osa_iter + timestep * inputScalars.subsets] % local_size[0];
			if (erotus[0] > 0)
				erotus[0] = (local_size[0] - erotus[0]);
			SET_RANGE2(localFP, local_size[0], 1);
			SET_LAUNCH_RANGE3(global, length[osa_iter + timestep * inputScalars.subsets] + erotus[0], 1, 1, localFP);
		}
#if defined(METAL)
		const uint32_t indD = osa_iter + timestep * inputScalars.subsets;
		kParams.d_N = { static_cast<unsigned int>(VEC_X(d_N[ii])), static_cast<unsigned int>(VEC_Y(d_N[ii])), static_cast<unsigned int>(VEC_Z(d_N[ii])) };
		kParams.d = d[ii];
		kParams.b = b[ii];
		kParams.d_bmax = bmax[ii];
		kParams.d_Scale4 = inputScalars.d_Scale4[ii];
		kParams.d_Scale5 = inputScalars.d_Scale[ii];
		kParams.dSize5 = inputScalars.dSize[ii];
		kParams.rings = inputScalars.rings;
		kParams.det_per_ring = inputScalars.det_per_ring;
		kParams.nProjections = length[indD];
		kParams.no_norm = no_norm;
		kParams.m_size = m_size;
		kParams.currentSubset = osa_iter;
		kParams.aa = ii;
		if (inputScalars.FPType == 2)
			kParams.orthWidth = inputScalars.tube_width;
		if (inputScalars.FPType == 3)
			kParams.orthWidth = inputScalars.cylRadiusProj3;
		SET_KERNEL_ARG_BYTES(encoder, kParams, sizeof(kParams), 0);
#endif // END METAL
		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", global[0]);
			mexPrintBase("global[1] = %u\n", global[1]);
			mexPrintBase("global[2] = %u\n", global[2]);
			mexPrintBase("localFP[0] = %u\n", localFP[0]);
			mexPrintBase("localFP[1] = %u\n", localFP[1]);
			mexPrintBase("localFP[2] = %u\n", localFP[2]);
			mexPrintBase("global[1] = %u\n", global[1]);
			mexPrintBase("global[2] = %u\n", global[2]);
			mexPrintBase("erotus[0] = %u\n", erotus[0]);
			mexPrintBase("erotus[1] = %u\n", erotus[1]);
			mexPrintBase("d_N[ii].s0 = %u\n", VEC_X(d_N[ii]));
			mexPrintBase("d_N[ii].s1 = %u\n", VEC_Y(d_N[ii]));
			mexPrintBase("d_N[ii].s2 = %u\n", VEC_Z(d_N[ii]));
			mexPrintBase("d[ii].s0 = %f\n", VEC_X(d[ii]));
			mexPrintBase("d[ii].s1 = %f\n", VEC_Y(d[ii]));
			mexPrintBase("d[ii].s2 = %f\n", VEC_Z(d[ii]));
			mexPrintBase("b[ii].s0 = %f\n", VEC_X(b[ii]));
			mexPrintBase("b[ii].s1 = %f\n", VEC_Y(b[ii]));
			mexPrintBase("b[ii].s2 = %f\n", VEC_Z(b[ii]));
			mexPrintBase("bmax[ii].s0 = %f\n", VEC_X(bmax[ii]));
			mexPrintBase("bmax[ii].s1 = %f\n", VEC_Y(bmax[ii]));
			mexPrintBase("bmax[ii].s2 = %f\n", VEC_Z(bmax[ii]));
#if defined(OPENCL)
			mexPrintBase("global.dimensions() = %u\n", global.dimensions());
			mexPrintBase("localFP.dimensions() = %u\n", localFP.dimensions());
#endif // END CUDA
			mexPrintBase("kernelIndFPSubIter = %u\n", kernelIndFPSubIter);
			mexPrintBase("kernelIndFP = %u\n", kernelIndFP);
			mexPrintBase("size_x = %u\n", inputScalars.nRowsD);
			mexPrintBase("size_y = %u\n", inputScalars.nColsD);
			mexPrintBase("m_size = %u\n", m_size);
			mexPrintBase("FPType = %u\n", inputScalars.FPType);
			if (inputScalars.FPType == 4) {
				mexPrintBase("dL = %f\n", inputScalars.dL);
				mexPrintBase("d_Scale4[ii].s[0] = %f\n", VEC_X(inputScalars.d_Scale4[ii]));
				mexPrintBase("d_Scale4[ii].s[1] = %f\n", VEC_Y(inputScalars.d_Scale4[ii]));
				mexPrintBase("d_Scale4[ii].s[2] = %f\n", VEC_Z(inputScalars.d_Scale4[ii]));
			}
			mexPrintBase("length[osa_iter] = %u\n", length[osa_iter + timestep * inputScalars.subsets]);
			mexPrintBase("listmode = %u\n", inputScalars.listmode);
			mexPrintBase("maskBP = %u\n", inputScalars.maskBP);
			mexPrintBase("maskFP = %u\n", inputScalars.maskFP);
			mexPrintBase("no_norm = %u\n", no_norm);
			mexPrintBase("ii = %u\n", ii);
			mexPrintBase("NVOXELS = %u\n", NVOXELS);
			mexPrintBase("NVOXELS5 = %u\n", NVOXELS5);
			mexPrintBase("osa_iter = %u\n", osa_iter);
			mexPrintBase("timestep = %u\n", timestep);
			mexPrintBase("memSize = %u\n", memSize);
			mexPrintBase("subsetType = %u\n", inputScalars.subsetType);
			mexEval();
		}
		TimerPoint tStart, tEnd;
		if (DEBUG || inputScalars.verbose >= 3) {
			FINISH_QUEUE(status, "", (STATUS_t)(-1));
			INIT_TIMER(tStart, tEnd);
		}
#if (defined(CUDA) || defined(HIP) || defined(OPENCL)) && !defined(AF)
		// ArrayFire and the projector share a stream; custom reconstruction inputs need an explicit synchronization.
		FINISH_QUEUE(status, "FP kernel failed\n", (STATUS_t)(-1));
#endif
#if defined(METAL)
		kernelIndFPSubIter = 1U;
		if (inputScalars.FPType >= 1 && inputScalars.FPType <= 3) {
			if (inputScalars.SPECT) {
				KARG_METAL_SLOT(kernelIndFPSubIter, 1);
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_rayShiftsDetector);
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_rayShiftsSource);
			}
			if (inputScalars.TOF) {
				KARG_METAL_SLOT(kernelIndFPSubIter, 3);
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_TOFCenter);
			}
			if (inputScalars.FPType == 2 || inputScalars.FPType == 3) {
				KARG_METAL_SLOT(kernelIndFPSubIter, 4);
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_V);
			}
		}
		else if (inputScalars.FPType == 4 && !inputScalars.CT && inputScalars.TOF) {
			KARG_METAL_SLOT(kernelIndFPSubIter, 1);
			KARG(kTemp, kernelFP, kernelIndFPSubIter, d_TOFCenter);
		}
#endif
		if (!inputScalars.CT && (inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3 || inputScalars.FPType == 4)) {
#if defined(METAL)
			if (inputScalars.FPType == 4)
				KARG_METAL_SLOT(kernelIndFPSubIter, 2);
			else
				KARG_METAL_SLOT(kernelIndFPSubIter, 5);
#endif
			if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation) {
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_atten[timestep][osa_iter]);
			}
			else if (inputScalars.attenuation_correction && inputScalars.CTAttenuation) {
				if (inputScalars.size_atten > inputScalars.im_dim[0]) {
					if (inputScalars.useBuffers) {
						KARG(kTemp, kernelFP, kernelIndFPSubIter, d_attenB[timestep]);
					}
					else
						KARG(kTemp, kernelFP, kernelIndFPSubIter, d_attenIm[timestep]);
				}
				else {
					if (inputScalars.useBuffers) {
						KARG(kTemp, kernelFP, kernelIndFPSubIter, d_attenB[0]);
					}
					else
						KARG(kTemp, kernelFP, kernelIndFPSubIter, d_attenIm[0]);
				}
			}
		}
		if (inputScalars.FPType == 5 || inputScalars.FPType == 4) {
			KARG_SCALAR(kTemp, kernelFP, kernelIndFPSubIter, d_N[ii]);
			KARG_SCALAR(kTemp, kernelFP, kernelIndFPSubIter, b[ii]);
			if (inputScalars.FPType == 5) {
				KARG_SCALAR(kTemp, kernelFP, kernelIndFPSubIter, inputScalars.dSize[ii]);
				KARG_SCALAR(kTemp, kernelFP, kernelIndFPSubIter, d[ii]);
				KARG_SCALAR(kTemp, kernelFP, kernelIndFPSubIter, inputScalars.d_Scale[ii]);
			}
			else {
				KARG_SCALAR(kTemp, kernelFP, kernelIndFPSubIter, bmax[ii]);
				KARG_SCALAR(kTemp, kernelFP, kernelIndFPSubIter, inputScalars.d_Scale4[ii]);
				if (inputScalars.largeDim) {
					// Only volume 0 is divided into subvolumes
					// for the multi-resolution volumes b/bmax already span the whole volume
					bzGlobalFP[0] = (ii == 0) ? inputScalars.lDimStruct.bz[0] : VEC_Z(b[ii]);
					bzGlobalFP[1] = (ii == 0) ? inputScalars.lDimStruct.bmaxZ[inputScalars.subsets - 1] : VEC_Z(bmax[ii]);
					KARG_SCALAR(kTemp, kernelFP, kernelIndFPSubIter, bzGlobalFP[0]);
					KARG_SCALAR(kTemp, kernelFP, kernelIndFPSubIter, bzGlobalFP[1]);
				}
			}
		}
		if (inputScalars.FPType == 4) {
			KARG_METAL_SLOT(kernelIndFPSubIter, 3);
#ifdef METAL
			if (updateMetalImageTextureFromBuffer(inputScalars, ii) != 0) {
				encoder->endEncoding();
				return -1;
			}
#endif
			KARG(kTemp, kernelFP, kernelIndFPSubIter, vec_opencl.d_image_os);
			KARG(kTemp, kernelFP, kernelIndFPSubIter, d_output);
			if (((inputScalars.listmode == 0 || inputScalars.indexBased) && !(inputScalars.CT || inputScalars.SPECT)) || 
				(!inputScalars.loadTOF && inputScalars.listmode > 0)) {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_x[0][0]);
			}
			else {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_x[timestep][osa_iter]);
			}
			if ((inputScalars.CT || inputScalars.PET)) {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_z[timestep][osa_iter]);
			}
			else if (inputScalars.listmode > 0) {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_z[0][0]);
			}
			else {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_z[timestep][inputScalars.osa_iter0]);
			}
			if (inputScalars.maskFP) {
				KARG_METAL_SLOT(kernelIndFPSubIter, 7);
				if (inputScalars.useBuffers) {
					int subset = 0;
					if (inputScalars.maskFPZ > 1 && !(inputScalars.SPECT && inputScalars.maskFPZ == inputScalars.nHeads))
						subset = osa_iter;
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_maskFPB[timestep][subset]);
				}
				else
					if (inputScalars.maskFPZ > 1) {
						KARG(kTemp, kernelFP, kernelIndFPSubIter, (inputScalars.SPECT && inputScalars.maskFPZ == inputScalars.nHeads) ? d_maskFP3[timestep][0] : d_maskFP3[timestep][osa_iter]);
					}
					else {
						KARG(kTemp, kernelFP, kernelIndFPSubIter, d_maskFP);
					}
			}
			KARG_SCALAR(kTemp, kernelFP, kernelIndFPSubIter, length[osa_iter + timestep * inputScalars.subsets]);
			if ((inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7) && inputScalars.subsets > 1 && !inputScalars.CT && !inputScalars.SPECT && !inputScalars.PET && inputScalars.listmode == 0) {
				KARG_METAL_SLOT(kernelIndFPSubIter, 9);
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_xyindex[osa_iter]);
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_zindex[osa_iter]);
			}
			if (inputScalars.listmode > 0 && inputScalars.indexBased) {
				KARG_METAL_SLOT(kernelIndFPSubIter, 9);
				if (!inputScalars.loadTOF) {
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_trIndex[0][0]);
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_axIndex[0][0]);
				}
				else {
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_trIndex[timestep][osa_iter]);
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_axIndex[timestep][osa_iter]);
				}
			}
			if (inputScalars.listmode > 0 && inputScalars.TOF) {
				KARG_METAL_SLOT(kernelIndFPSubIter, 11);
				if (!inputScalars.loadTOF) {
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_TOFIndex[0][0]);
				}
				else {
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_TOFIndex[timestep][osa_iter]);
				}
			}
			if (inputScalars.raw) {
				KARG_METAL_SLOT(kernelIndFPSubIter, 12);
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_L[osa_iter]);
				KARG_SCALAR(kTemp, kernelFP, kernelIndFPSubIter, inputScalars.det_per_ring);
			}
			if (inputScalars.normalization_correction) {
				KARG_METAL_SLOT(kernelIndFPSubIter, 13);
				// TODO: Listmode normalization
				//if (inputScalars.listmode > 0 && inputScalars.indexBased)
				//	status = kernelFP.setArg(kernelIndFPSubIter++, d_norm[0]);
				//else
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_norm[timestep][osa_iter]);
			}
			if (inputScalars.scatter) {
				KARG_METAL_SLOT(kernelIndFPSubIter, 14);
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_scat[timestep][osa_iter]);
			}
			KARG_SCALAR(kTemp, kernelFP, kernelIndFPSubIter, no_norm);
			KARG_SCALAR(kTemp, kernelFP, kernelIndFPSubIter, m_size);
			KARG_SCALAR(kTemp, kernelFP, kernelIndFPSubIter, osa_iter);
			KARG_SCALAR(kTemp, kernelFP, kernelIndFPSubIter, ii);
		}
		else if (inputScalars.FPType == 5) {
			KARG_METAL_SLOT(kernelIndFPSubIter, 1);
			if (!inputScalars.loadTOF && inputScalars.listmode > 0) {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_x[0][0]);
			}
			else
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_x[timestep][osa_iter]);
			KARG(kTemp, kernelFP, kernelIndFPSubIter, d_z[timestep][osa_iter]);
			KARG(kTemp, kernelFP, kernelIndFPSubIter, vec_opencl.d_image_os);
			KARG(kTemp, kernelFP, kernelIndFPSubIter, vec_opencl.d_image_os_int);
			KARG(kTemp, kernelFP, kernelIndFPSubIter, d_output);
			if (inputScalars.meanFP) {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_meanFP);
			}
			if (inputScalars.maskFP) {
				KARG_METAL_SLOT(kernelIndFPSubIter, 7);
				if (inputScalars.useBuffers) {
					int subset = 0;
					if (inputScalars.maskFPZ > 1 && !(inputScalars.SPECT && inputScalars.maskFPZ == inputScalars.nHeads))
						subset = osa_iter;
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_maskFPB[timestep][subset]);
				}
				else
					if (inputScalars.maskFPZ > 1) {
						KARG(kTemp, kernelFP, kernelIndFPSubIter, (inputScalars.SPECT && inputScalars.maskFPZ == inputScalars.nHeads) ? d_maskFP3[timestep][0] : d_maskFP3[timestep][osa_iter]);
					}
					else
						KARG(kTemp, kernelFP, kernelIndFPSubIter, d_maskFP);
			}
			if (inputScalars.normalization_correction) {
				// TODO: Listmode normalization
				//if (inputScalars.listmode > 0 && inputScalars.indexBased)
				//	status = kernelFP.setArg(kernelIndFPSubIter++, d_norm[0]);
				//else
				KARG_METAL_SLOT(kernelIndFPSubIter, 8);
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_norm[timestep][osa_iter]);
			}
			KARG_SCALAR(kTemp, kernelFP, kernelIndFPSubIter, length[osa_iter + timestep * inputScalars.subsets]);
		}
		else if ((inputScalars.FPType == 1 || inputScalars.FPType == 2 || inputScalars.FPType == 3)) {
			if (inputScalars.maskFP) {
				KARG_METAL_SLOT(kernelIndFPSubIter, 6);
				if (inputScalars.useBuffers) {
					int subset = 0;
					if (inputScalars.maskFPZ > 1 && !(inputScalars.SPECT && inputScalars.maskFPZ == inputScalars.nHeads))
						subset = osa_iter;
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_maskFPB[timestep][subset]);
				}
				else
					if (inputScalars.maskFPZ > 1) {
						KARG(kTemp, kernelFP, kernelIndFPSubIter, (inputScalars.SPECT && inputScalars.maskFPZ == inputScalars.nHeads) ? d_maskFP3[timestep][0] : d_maskFP3[timestep][osa_iter]);
					}
					else
						KARG(kTemp, kernelFP, kernelIndFPSubIter, d_maskFP);
			}
			if ((inputScalars.CT || inputScalars.PET || inputScalars.SPECT) && inputScalars.listmode == 0) {
				KARG_SCALAR(kTemp, kernelFP, kernelIndFPSubIter, length[osa_iter + timestep * inputScalars.subsets]);
			}
			KARG_METAL_SLOT(kernelIndFPSubIter, 8);
			if (((inputScalars.listmode == 0 || inputScalars.indexBased) && !(inputScalars.CT || inputScalars.SPECT)) || (!inputScalars.loadTOF && inputScalars.listmode > 0)) {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_x[0][0]);
			}
			else
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_x[timestep][osa_iter]);
			if ((inputScalars.CT || inputScalars.PET || inputScalars.SPECT)) {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_z[timestep][osa_iter]);
			}
			else if (inputScalars.listmode > 0) {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_z[0][0]);
			}
			else
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_z[timestep][inputScalars.osa_iter0]);
			if (inputScalars.normalization_correction) {
				// TODO: Listmode normalization
				//if (inputScalars.listmode > 0 && inputScalars.indexBased)
				//	status = kernelFP.setArg(kernelIndFPSubIter++, d_norm[0]);
				//else
				KARG_METAL_SLOT(kernelIndFPSubIter, 10);
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_norm[timestep][osa_iter]);
			}
			if (inputScalars.scatter) {
				KARG_METAL_SLOT(kernelIndFPSubIter, 11);
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_scat[timestep][osa_iter]);
			}
			KARG_METAL_SLOT(kernelIndFPSubIter, 12);
			KARG(kTemp, kernelFP, kernelIndFPSubIter, d_Summ[uu]);
			KARG_SCALAR(kTemp, kernelFP, kernelIndFPSubIter, d_N[ii]);
			KARG_SCALAR(kTemp, kernelFP, kernelIndFPSubIter, d[ii]);
			KARG_SCALAR(kTemp, kernelFP, kernelIndFPSubIter, b[ii]);
			KARG_SCALAR(kTemp, kernelFP, kernelIndFPSubIter, bmax[ii]);
			if ((inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7) && inputScalars.subsets > 1 && !inputScalars.CT && !inputScalars.SPECT && !inputScalars.PET && inputScalars.listmode == 0) {
				KARG_METAL_SLOT(kernelIndFPSubIter, 13);
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_xyindex[osa_iter]);
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_zindex[osa_iter]);
			}
			if (inputScalars.listmode > 0 && inputScalars.indexBased) {
				KARG_METAL_SLOT(kernelIndFPSubIter, 15);
				if (!inputScalars.loadTOF) {
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_trIndex[0][0]);
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_axIndex[0][0]);
				}
				else {
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_trIndex[timestep][osa_iter]);
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_axIndex[timestep][osa_iter]);
				}
			}
			if (inputScalars.listmode > 0 && inputScalars.TOF) {
				KARG_METAL_SLOT(kernelIndFPSubIter, 17);
				if (!inputScalars.loadTOF) {
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_TOFIndex[0][0]);
				}
				else {
					KARG(kTemp, kernelFP, kernelIndFPSubIter, d_TOFIndex[timestep][osa_iter]);
				}
			}
			if (inputScalars.raw) {
				KARG_METAL_SLOT(kernelIndFPSubIter, 18);
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_L[osa_iter]);
			}
			KARG_METAL_SLOT(kernelIndFPSubIter, 19);
			if (inputScalars.useBuffers) {
				KARG(kTemp, kernelFP, kernelIndFPSubIter, vec_opencl.d_im);
			}
				else {
#if defined(METAL)
					// The ArrayFire image estimate changes after every subset and
					// d_im also changes between multiresolution volumes. Refresh
					// the Metal texture for every forward projection.
					if (updateMetalImageTextureFromBuffer(inputScalars, ii) != 0) {
						encoder->endEncoding();
						return -1;
					}
#endif
					KARG(kTemp, kernelFP, kernelIndFPSubIter, vec_opencl.d_image_os);
				}
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_output);
			if (inputScalars.SPECT) {
				KARG_METAL_SLOT(kernelIndFPSubIter, 21);
				KARG(kTemp, kernelFP, kernelIndFPSubIter, d_detectorVector[timestep][osa_iter]);
			}
			KARG_SCALAR(kTemp, kernelFP, kernelIndFPSubIter, no_norm);
			KARG_SCALAR(kTemp, kernelFP, kernelIndFPSubIter, m_size);
			KARG_SCALAR(kTemp, kernelFP, kernelIndFPSubIter, osa_iter);
			KARG_SCALAR(kTemp, kernelFP, kernelIndFPSubIter, ii);
		}
		if (DEBUG || inputScalars.verbose >= 3)
			START_TIMER(tStart);
#if defined(CUDA) || defined(HIP)
		status = cuLaunchKernel(kernelFP, global[0], global[1], global[2], localFP[0], localFP[1], localFP[2], 0, CLCommandQueue[0], kTemp.data(), NULL);
		CUDA_CHECK(status, "Failed to launch forward projection kernel\n", -1);
#elif defined(METAL)
		{
			const MTL::Size threadsPerThreadgroup = MTL::Size::Make(localFP[0], localFP[1], localFP[2]);
			const MTL::Size threadgroupsPerGrid = MTL::Size::Make(global[0] / localFP[0], global[1] / localFP[1], global[2] / localFP[2]);
			encoder->dispatchThreadgroups(threadgroupsPerGrid, threadsPerThreadgroup);
			encoder->endEncoding();
			commandBuffer->commit();
			commandBuffer->waitUntilCompleted();
		}
#elif defined(OPENCL)
		status = CLCommandQueue[0].enqueueNDRangeKernel(kernelFP, cl::NDRange(), global, localFP, NULL);
		OCL_CHECK(status, "\n", -1);
#endif // END CUDA
		if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Forward projection kernel launched successfully\n");
		}
#if (defined(CUDA) || defined(HIP) || defined(OPENCL)) && !defined(AF)
		FINISH_QUEUE(status, "FP kernel failed to complete\n", (STATUS_t)(-1));
#endif
		if (DEBUG || inputScalars.verbose >= 3) {
#if defined(OPENCL)
			CLCommandQueue[0].finish();
			STOP_TIMER(tEnd);
#elif defined(METAL) || defined(CUDA) || defined(HIP)
			STOP_TIMER(tEnd);
#endif
			PRINT_TIMER(tStart, tEnd, "Forward projection completed in %f seconds\n");
#if defined(CUDA) || defined(HIP)
			cuEventDestroy(tStart);
			cuEventDestroy(tEnd);
#endif
		}
		return 0;
	}

	/// <summary>
	/// Compute the backprojection for the selected projector type
	/// </summary>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="w_vec specifies some of the special options/parameters used"></param>
	/// <param name="osa_iter current subset (sub-iteration)"></param>
	/// <param name="length the number of measurements/projection/sinograms per subset"></param>
	/// <param name="m_size for projector types 1-3, the total number of LORs"></param>
	/// <param name="compSens if true, computes the sensitivity image as well"></param>
	/// <returns></returns>
#if defined(CUDA) || defined(HIP)
	inline int backwardProjection(scalarStruct & inputScalars, Weighting & w_vec, uint32_t osa_iter, uint32_t timestep,
		std::vector<int64_t>&length, uint64_t m_size, const RecMethods & MethodList = RecMethods(), const bool compSens = false, int ii = 0, const int uu = 0,
		int ee = -1, const int queueIdx = 0, const bool newInput = true) {
#elif defined(METAL) || defined(OPENCL)
	// MethodList defaults to an empty RecMethods() so the METAL branch below keeps compiling unchanged
	// TODO: Metal support for fastPDHG?
	inline int backwardProjection(const scalarStruct & inputScalars, Weighting & w_vec, const uint32_t osa_iter, const uint32_t timestep,
		const std::vector<int64_t>&length, const uint64_t m_size, const RecMethods & MethodList = RecMethods(), const bool compSens = false, const int32_t ii = 0, const int uu = 0,
		int ee = -1, const int queueIdx = 0, const bool newInput = true) {
#endif // END CUDA
		if (inputScalars.verbose >= 3 || DEBUG)
		    mexPrintVar("Starting backprojection for projector type = ", inputScalars.BPType);
		STATUS_t status = SUCCESS_VALUE;
#if defined(CUDA) || defined(HIP)
		std::vector<void*> kTemp = BPArgs;
#elif defined(OPENCL)
		kernelIndBPSubIter = kernelIndBP;
#endif // END CUDA
        if (ee < 0)
			ee = uu;
		if (inputScalars.listmode > 0 && compSens) {
			kernelApu = kernelBP;
			kernelBP = kernelSensList;
#if defined(CUDA) || defined(HIP)
			kTemp = SensArgs;
#elif defined(OPENCL)
			kernelIndBPSubIter = kernelIndSens;
#endif // END CUDA
		}
		int indD = osa_iter + timestep * inputScalars.subsets;
#if defined(METAL)
		const bool useMetalSideQueue = queueIdx > 0 && static_cast<size_t>(queueIdx) <= sideQueues.size();
		MTL::CommandQueue* activeQueue = useMetalSideQueue ? sideQueues[queueIdx - 1].get() : queueBP.get();
		if (!activeQueue || !kernelBP) {
			mexPrint("Unable to create Metal backprojection encoder");
			return -1;
		}
		NS::SharedPtr<MTL::CommandBuffer> commandBuffer = NS::RetainPtr(activeQueue->commandBuffer());
		NS::SharedPtr<MTL::ComputeCommandEncoder> encoder = NS::RetainPtr(commandBuffer->computeCommandEncoder());
		if (!commandBuffer || !encoder) {
			mexPrint("Unable to create Metal backprojection encoder");
			return -1;
		}
		encoder->setComputePipelineState(kernelBP.get());
#endif // END METAL
		WORKRANGE_t localBP = local;
		TimerPoint tStart, tEnd;
		if (DEBUG || inputScalars.verbose >= 3) {
			INIT_TIMER(tStart, tEnd);
		}

#if defined(METAL)
		kernelIndBPSubIter = 1U;
		if (inputScalars.BPType >= 1 && inputScalars.BPType <= 3) {
			if (inputScalars.SPECT) {
				KARG_METAL_SLOT(kernelIndBPSubIter, 1);
				KARG(kTemp, kernelBP, kernelIndBPSubIter, d_rayShiftsDetector);
				KARG(kTemp, kernelBP, kernelIndBPSubIter, d_rayShiftsSource);
			}
			if (inputScalars.TOF) {
				KARG_METAL_SLOT(kernelIndBPSubIter, 3);
				KARG(kTemp, kernelBP, kernelIndBPSubIter, d_TOFCenter);
			}
			if (inputScalars.BPType == 2 || inputScalars.BPType == 3) {
				KARG_METAL_SLOT(kernelIndBPSubIter, 4);
				KARG(kTemp, kernelBP, kernelIndBPSubIter, d_V);
			}
		}
		else if (inputScalars.BPType == 4 && !inputScalars.CT && inputScalars.TOF) {
			KARG_METAL_SLOT(kernelIndBPSubIter, 1);
			KARG(kTemp, kernelBP, kernelIndBPSubIter, d_TOFCenter);
		}
		kParams.d_N = { static_cast<unsigned int>(VEC_X(d_N[ii])), static_cast<unsigned int>(VEC_Y(d_N[ii])), static_cast<unsigned int>(VEC_Z(d_N[ii])) };
		kParams.d = d[ii];
		kParams.b = b[ii];
		kParams.d_bmax = bmax[ii];
		kParams.d_Scale4 = inputScalars.d_Scale4[ii];
		kParams.d_Scale5 = inputScalars.d_Scale[ii];
		kParams.dSize5 = inputScalars.dSizeBP;
		kParams.kerroin4 = (inputScalars.BPType == 4 && w_vec.kerroin4) ? w_vec.kerroin4[ii] : 0.f;
		kParams.nProjections = length[indD];
		kParams.no_norm = no_norm;
		kParams.m_size = m_size;
		kParams.currentSubset = osa_iter;
		kParams.aa = ii;
		if (inputScalars.BPType == 2)
			kParams.orthWidth = inputScalars.tube_width;
		if (inputScalars.BPType == 3)
			kParams.orthWidth = inputScalars.cylRadiusProj3;
		SET_KERNEL_ARG_BYTES(encoder, kParams, sizeof(kParams), 0);
#endif
		if (!inputScalars.CT && (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3 || inputScalars.BPType == 4)) {
#if defined(METAL)
			if (inputScalars.BPType == 4)
				KARG_METAL_SLOT(kernelIndBPSubIter, 2);
			else
				KARG_METAL_SLOT(kernelIndBPSubIter, 5);
#endif
			if (inputScalars.attenuation_correction && !inputScalars.CTAttenuation) {
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_atten[timestep][osa_iter]);
			}
			else if (inputScalars.attenuation_correction && inputScalars.CTAttenuation) {
				if (inputScalars.size_atten > inputScalars.im_dim[0]) {
					if (inputScalars.useBuffers) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_attenB[timestep]);
					}
					else
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_attenIm[timestep]);
				}
				else {
					if (inputScalars.useBuffers) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_attenB[0]);
					}
					else
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_attenIm[0]);
				}
			}
		}
		if (inputScalars.BPType == 1 || inputScalars.BPType == 2 || inputScalars.BPType == 3) {
			if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0) {
				SET_LAUNCH_RANGE3(global, inputScalars.nRowsD + erotus[0], inputScalars.nColsD + erotus[1], length[indD], local);
			}
			else if (inputScalars.listmode > 0 && compSens) {
				const size_t sensitivityDepth = static_cast<size_t>(inputScalars.rings) * static_cast<size_t>(inputScalars.rings) * static_cast<size_t>(inputScalars.nLayers);
				SET_LAUNCH_RANGE3(global, static_cast<size_t>(inputScalars.det_per_ring) + erotusSens[0], 
					static_cast<size_t>(inputScalars.det_per_ring) + erotusSens[1], sensitivityDepth, local);
			}
			else {
				erotus[0] = length[indD] % local_size[0];
				if (erotus[0] > 0)
					erotus[0] = (local_size[0] - erotus[0]);
				SET_LAUNCH_RANGE3(global, length[indD] + erotus[0], 1, 1, local);
			}
			if (DEBUG) {
				mexPrintBase("global[0] = %u\n", global[0]);
				mexPrintBase("localBP[0] = %u\n", localBP[0]);
				mexPrintBase("localBP[1] = %u\n", localBP[1]);
				mexPrintBase("global[1] = %u\n", global[1]);
				mexPrintBase("global[2] = %u\n", global[2]);
				if (inputScalars.listmode > 0 && compSens) {
					mexPrintBase("erotusSens[0] = %u\n", erotusSens[0]);
					mexPrintBase("erotusSens[1] = %u\n", erotusSens[1]);
				}
				else {
					mexPrintBase("erotus[0] = %u\n", erotus[0]);
					mexPrintBase("erotus[1] = %u\n", erotus[1]);
				}
#if defined(CUDA) || defined(HIP)
				mexPrintBase("kernelIndBPSubIter = %u\n", BPArgs.size());
#elif defined(OPENCL)
				mexPrintBase("global.dimensions() = %u\n", global.dimensions());
				mexPrintBase("localBP.dimensions() = %u\n", localBP.dimensions());
				mexPrintBase("kernelIndBPSubIter = %u\n", kernelIndBPSubIter);
#endif // END CUDA
				mexPrintBase("m_size = %u\n", m_size);
				mexPrintBase("nRowsD = %u\n", inputScalars.nRowsD);
				mexPrintBase("nColsD = %u\n", inputScalars.nColsD);
				mexPrintBase("length[indD] = %u\n", length[indD]);
				mexPrintBase("listmode = %u\n", inputScalars.listmode);
				mexPrintBase("rings = %u\n", inputScalars.rings);
				mexPrintBase("im_dim = %u\n", inputScalars.im_dim[ii]);
				mexPrintBase("no_norm = %u\n", no_norm);
				mexPrintBase("ii = %u\n", ii);
				mexPrintBase("memSize = %u\n", memSize);
				mexPrintBase("compSens = %u\n", compSens);
				mexPrintBase("osa_iter = %u\n", osa_iter);
				mexPrintBase("timestep = %u\n", timestep);
				mexEval();
			}
			// Set kernelBP arguments
			if (inputScalars.maskFP || inputScalars.maskBP) {
				if (inputScalars.maskFP) {
					KARG_METAL_SLOT(kernelIndBPSubIter, 6);
					if (inputScalars.useBuffers) {
						int subset = 0;
						if (inputScalars.maskFPZ > 1 && !(inputScalars.SPECT && inputScalars.maskFPZ == inputScalars.nHeads))
							subset = osa_iter;
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskFPB[timestep][subset]);
					}
					else
						if (inputScalars.maskFPZ > 1) {
							KARG(kTemp, kernelBP, kernelIndBPSubIter, (inputScalars.SPECT && inputScalars.maskFPZ == inputScalars.nHeads) ? d_maskFP3[timestep][0] : d_maskFP3[timestep][osa_iter]);
						}
						else
							KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskFP);
				}
				if (inputScalars.maskBP) {
					KARG_METAL_SLOT(kernelIndBPSubIter, 7);
					if (inputScalars.useBuffers) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskBPB[timestep][ii]);
					}
					else {
						if (inputScalars.maskBPZ > 1) {
							KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskBP3[timestep][ii]);
						}
						else
							KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskBP[timestep][ii]);
					}
				}
			}
			if ((inputScalars.CT || inputScalars.PET || inputScalars.SPECT) && inputScalars.listmode == 0)
				KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, length[indD]);
			KARG_METAL_SLOT(kernelIndBPSubIter, 8);
			if (compSens) {
				KARG(kTemp, kernelBP, kernelIndBPSubIter, d_xFull[0]);
				KARG(kTemp, kernelBP, kernelIndBPSubIter, d_zFull[0]);
				KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, inputScalars.rings);
			}
			else {
				if (((inputScalars.listmode == 0 || inputScalars.indexBased) && !(inputScalars.CT || inputScalars.SPECT)) || (!inputScalars.loadTOF && inputScalars.listmode > 0)) {
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_x[0][0]);
				}
				else
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_x[timestep][osa_iter]);
				if ((inputScalars.CT || inputScalars.PET || inputScalars.SPECT)) {
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_z[timestep][osa_iter]);
				}
				else if (inputScalars.listmode > 0) {
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_z[0][0]);
				}
				else
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_z[timestep][inputScalars.osa_iter0]);
			}
			if (compSens) {
				if (inputScalars.normalization_correction) {
					KARG_METAL_SLOT(kernelIndBPSubIter, 10);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_normFull[0]);
				}
				if (inputScalars.scatter) {
					KARG_METAL_SLOT(kernelIndBPSubIter, 11);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_scatFull[0]);
				}
			}
			else {
				if (inputScalars.normalization_correction) {
					KARG_METAL_SLOT(kernelIndBPSubIter, 10);
					// TODO: Listmode normalization
					//if (inputScalars.listmode > 0 && inputScalars.indexBased)
					//	status = kernelBP.setArg(kernelIndBPSubIter++, d_norm[0]);
					//else
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_norm[timestep][osa_iter]);
				}
				if (inputScalars.scatter) {
					KARG_METAL_SLOT(kernelIndBPSubIter, 11);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_scat[timestep][osa_iter]);
				}
			}
			KARG_METAL_SLOT(kernelIndBPSubIter, 12);
			KARG(kTemp, kernelBP, kernelIndBPSubIter, d_Summ[ee]);
			KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, d_N[ii]);
			KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, d[ii]);
			KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, b[ii]);
			KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, bmax[ii]);
			if ((inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7) && inputScalars.subsets > 1 && !inputScalars.CT && !inputScalars.SPECT && !inputScalars.PET && inputScalars.listmode == 0) {
				KARG_METAL_SLOT(kernelIndBPSubIter, 13);
				KARG(kTemp, kernelBP, kernelIndBPSubIter, d_xyindex[osa_iter]);
				KARG(kTemp, kernelBP, kernelIndBPSubIter, d_zindex[osa_iter]);
			}
			if (inputScalars.listmode > 0 && inputScalars.indexBased && !compSens) {
				KARG_METAL_SLOT(kernelIndBPSubIter, 15);
				if (!inputScalars.loadTOF) {
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_trIndex[0][0]);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_axIndex[0][0]);
				}
				else {
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_trIndex[timestep][osa_iter]);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_axIndex[timestep][osa_iter]);
				}
			}
			if (inputScalars.listmode > 0 && inputScalars.TOF) {
				KARG_METAL_SLOT(kernelIndBPSubIter, 17);
				if (!inputScalars.loadTOF) {
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_TOFIndex[0][0]);
				}
				else {
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_TOFIndex[timestep][osa_iter]);
				}
			}
			// TODO: Raw data?
			if (inputScalars.raw) {
				KARG_METAL_SLOT(kernelIndBPSubIter, 18);
				KARG(kTemp, kernelBP, kernelIndBPSubIter, d_L[osa_iter]);
			}
			KARG_METAL_SLOT(kernelIndBPSubIter, 19);
			KARG(kTemp, kernelBP, kernelIndBPSubIter, d_output);
			KARG(kTemp, kernelBP, kernelIndBPSubIter, vec_opencl.d_rhs_os[uu]);
			if (inputScalars.SPECT) {
				KARG_METAL_SLOT(kernelIndBPSubIter, 21);
				KARG(kTemp, kernelBP, kernelIndBPSubIter, d_detectorVector[timestep][osa_iter]);
			}
			KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, no_norm);
			KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, m_size);
			KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, osa_iter);
			KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, ii);
			}
		else {
			if (inputScalars.CT) {
				// Projector type 7 always reads the projection data from a buffer, no image/texture is needed
				if (!inputScalars.useBuffers && inputScalars.BPType != 7) {
#if defined(CUDA) || defined(HIP)
					CUDA_ARRAY3D_DESCRIPTOR_st arr3DDesc;
					std::memset(&arr3DDesc, 0, sizeof(arr3DDesc));
					arr3DDesc.Format = CUarray_format::CU_AD_FORMAT_FLOAT;
					arr3DDesc.NumChannels = 1;
					arr3DDesc.Height = inputScalars.nColsD;
					arr3DDesc.Width = inputScalars.nRowsD;
					arr3DDesc.Depth = length[indD];
					if (inputScalars.BPType == 5) {
						arr3DDesc.Height++;
						arr3DDesc.Width++;
					}
					if (DEBUG) {
						mexPrintBase("arr3DDesc.Height = %u\n", arr3DDesc.Height);
						mexPrintBase("arr3DDesc.Width = %u\n", arr3DDesc.Width);
						mexPrintBase("arr3DDesc.Depth = %u\n", arr3DDesc.Depth);
						mexEval();
					}
					if (newInput) {
						if (BPImageDims[0] != arr3DDesc.Width || BPImageDims[1] != arr3DDesc.Height || BPImageDims[2] != arr3DDesc.Depth) {
							if (BPImageDims[2] != 0) {
								getErrorString(cuTexObjectDestroy(d_inputImage));
								getErrorString(cuArrayDestroy(BPArray));
							}
							status = cuArray3DCreate(&BPArray, &arr3DDesc);
							CUDA_CHECK(status, "Array creation failed\n", -1);
							CUDA_TEXTURE_DESC texDesc;
							CUDA_RESOURCE_DESC resDescIm;
							CUDA_RESOURCE_VIEW_DESC viewDesc;
							std::memset(&texDesc, 0, sizeof(texDesc));
							std::memset(&resDescIm, 0, sizeof(resDescIm));
							std::memset(&viewDesc, 0, sizeof(viewDesc));
							viewDesc.height = inputScalars.nColsD;
							viewDesc.width = inputScalars.nRowsD;
							viewDesc.depth = length[indD];
							viewDesc.format = CUresourceViewFormat::CU_RES_VIEW_FORMAT_FLOAT_1X32;
							resDescIm.resType = CUresourcetype::CU_RESOURCE_TYPE_ARRAY;
							resDescIm.res.array.hArray = BPArray;
							texDesc.addressMode[0] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
							texDesc.addressMode[1] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
							texDesc.addressMode[2] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
							texDesc.flags = CU_TRSF_NORMALIZED_COORDINATES;
							texDesc.filterMode = CUfilter_mode::CU_TR_FILTER_MODE_LINEAR;
							if (inputScalars.BPType == 5) {
								viewDesc.height++;
								viewDesc.width++;
							}
							status = cuTexObjectCreate(&d_inputImage, &resDescIm, &texDesc, &viewDesc);
							CUDA_CHECK(status, "Image creation failed\n", -1);
							BPImageDims[0] = arr3DDesc.Width;
							BPImageDims[1] = arr3DDesc.Height;
							BPImageDims[2] = arr3DDesc.Depth;
						}
						CUDA_MEMCPY3D cpy3d;
						std::memset(&cpy3d, 0, sizeof(cpy3d));
						cpy3d.srcMemoryType = CUmemorytype::CU_MEMORYTYPE_DEVICE;
						cpy3d.srcDevice = reinterpret_cast<CUdeviceptr>(d_output);
						cpy3d.srcPitch = inputScalars.nRowsD * sizeof(float);
						cpy3d.srcHeight = inputScalars.nColsD;
						cpy3d.dstMemoryType = CUmemorytype::CU_MEMORYTYPE_ARRAY;
						cpy3d.dstArray = BPArray;
						cpy3d.WidthInBytes = inputScalars.nRowsD * sizeof(float);
						cpy3d.Height = inputScalars.nColsD;
						cpy3d.Depth = length[indD];
						if (inputScalars.BPType == 5) {
							cpy3d.srcPitch += sizeof(float);
							cpy3d.srcHeight++;
							cpy3d.WidthInBytes += sizeof(float);
							cpy3d.Height++;
						}
						status = cuMemcpy3DAsync(&cpy3d, CLCommandQueue[0]);
						CUDA_CHECK(status, "Array mem copy failed\n", -1);
					}
#elif defined(OPENCL)
					cl::size_type imX = inputScalars.nRowsD;
					cl::size_type imY = inputScalars.nColsD;
					cl::size_type imZ = length[indD];
					if (inputScalars.BPType == 5) {
						imX++;
						imY++;
					}
					if (DEBUG) {
						mexPrintBase("image width = %u\n", imX);
						mexPrintBase("image height = %u\n", imY);
						mexPrintBase("image depth = %u\n", imZ);
						mexEval();
					}
					cl::detail::size_t_array region = { imX, imY, imZ };
					// Make the inputs persistent such that they are only created when required
					if (newInput) {
						if (BPImageDims[0] != imX || BPImageDims[1] != imY || BPImageDims[2] != imZ) {
							d_inputImage = cl::Image3D(CLContext, CL_MEM_READ_ONLY, format, imX, imY, imZ, 0, 0, NULL, &status);
							OCL_CHECK(status, "Image creation failed\n", -1);
							BPImageDims[0] = imX;
							BPImageDims[1] = imY;
							BPImageDims[2] = imZ;
						}
						status = CLCommandQueue[0].enqueueCopyBufferToImage(d_output, d_inputImage, 0, origin, region);
						OCL_CHECK(status, "Image copy failed\n", -1);
					}
				#elif defined(METAL)
					size_t textureHeight = inputScalars.nColsD;
					size_t textureWidth = inputScalars.nRowsD;
					const size_t textureDepth = length[indD];
					if (inputScalars.BPType == 5) {
						textureHeight++;
						textureWidth++;
					}
					if (newInput) {
						if (!d_output || !d_output->contents()) {
							mexPrint("Metal backprojection input buffer is unavailable\n");
							return -1;
						}
						if (!d_inputImage || BPImageDims[0] != textureWidth || BPImageDims[1] != textureHeight || BPImageDims[2] != textureDepth) {
							d_inputImage = createMetalFloatTextureEmpty(
								metalTextureSpec(textureWidth, textureHeight, textureDepth, true));
							if (!d_inputImage) {
								mexPrint("Metal backprojection input texture creation failed\n");
								return -1;
							}
							BPImageDims[0] = textureWidth;
							BPImageDims[1] = textureHeight;
							BPImageDims[2] = textureDepth;
						}
						const MTL::Region textureRegion(0, 0, 0,
							static_cast<NS::UInteger>(textureWidth),
							static_cast<NS::UInteger>(textureHeight),
							static_cast<NS::UInteger>(textureDepth));
						const NS::UInteger bytesPerRow = static_cast<NS::UInteger>(textureWidth * sizeof(float));
						const NS::UInteger bytesPerImage = bytesPerRow * static_cast<NS::UInteger>(textureHeight);
						d_inputImage->replaceRegion(textureRegion, 0, 0, d_output->contents(), bytesPerRow, bytesPerImage);
					}
					else if (!d_inputImage) {
						mexPrint("Metal backprojection input texture cannot be reused before it is initialized\n");
						return -1;
					}
#endif // END CUDA
				}
#if defined(CUDA) || defined(HIP)
			if (inputScalars.BPType == 4) {
				global[0] = (inputScalars.Nx[ii] + erotusBP[0][ii]) / localBP[0];
				global[1] = (inputScalars.Ny[ii] + erotusBP[1][ii]) / localBP[1];
#elif defined(OPENCL) || defined(METAL)
				if (inputScalars.BPType == 4)
#endif // END CUDA
					if (!inputScalars.largeDim)
						if (!inputScalars.useHelical)
#if defined(CUDA) || defined(HIP)
							global[2] = (inputScalars.Nz[ii] + NVOXELS - 1) / NVOXELS;
#elif defined(OPENCL) || defined(METAL)
							SET_RANGE3(global, inputScalars.Nx[ii] + erotusBP[0][ii], inputScalars.Ny[ii] + erotusBP[1][ii], (inputScalars.Nz[ii] + NVOXELS - 1) / NVOXELS);
#endif // END CUDA
						else
#if defined(CUDA) || defined(HIP)
							global[2] = (inputScalars.Nz[ii] + NVOXELSHELICAL - 1) / NVOXELSHELICAL;
#elif defined(OPENCL) || defined(METAL)
							SET_RANGE3(global, inputScalars.Nx[ii] + erotusBP[0][ii], inputScalars.Ny[ii] + erotusBP[1][ii], (inputScalars.Nz[ii] + NVOXELSHELICAL - 1) / NVOXELSHELICAL);
#endif // END CUDA
					else
#if defined(CUDA) || defined(HIP)
						global[2] = inputScalars.Nz[ii];
			}
#elif defined(OPENCL) || defined(METAL)
						SET_RANGE3(global, inputScalars.Nx[ii] + erotusBP[0][ii], inputScalars.Ny[ii] + erotusBP[1][ii], inputScalars.Nz[ii]);
#endif // END CUDA
				else if (inputScalars.BPType == 5) {
#if defined(CUDA) || defined(HIP)
					if (inputScalars.pitch) {
						global[0] = (inputScalars.Nx[ii] + erotusBP[0][ii]) / localBP[0];
						global[1] = (inputScalars.Ny[ii] + erotusBP[1][ii]) / localBP[1];
						global[2] = inputScalars.Nz[ii];
#elif defined(OPENCL) || defined(METAL)
					if (inputScalars.pitch)
						SET_RANGE3(global, inputScalars.Nx[ii] + erotusBP[0][ii], inputScalars.Ny[ii] + erotusBP[1][ii], inputScalars.Nz[ii]);
					else
						SET_RANGE3(global, inputScalars.Nx[ii] + erotusBP[0][ii], inputScalars.Ny[ii] + erotusBP[1][ii], (inputScalars.Nz[ii] + NVOXELS5 - 1) / NVOXELS5);
#endif // END CUDA
					}
#if defined(CUDA) || defined(HIP)
				else {
					global[0] = (inputScalars.Nx[ii] + erotusBP[0][ii]) / localBP[0];
					global[1] = (inputScalars.Ny[ii] + erotusBP[1][ii]) / localBP[1];
					global[2] = (inputScalars.Nz[ii] + NVOXELS5 - 1) / NVOXELS5;
				}
				}
			else {
				global[0] = (inputScalars.Nx[ii] + erotusBP[0][ii]) / localBP[0];
				global[1] = (inputScalars.Ny[ii] + erotusBP[1][ii]) / localBP[1];
				global[2] = inputScalars.Nz[ii];
			}
#elif defined(OPENCL) || defined(METAL)
				else
					SET_RANGE3(global, inputScalars.Nx[ii] + erotusBP[0][ii], inputScalars.Ny[ii] + erotusBP[1][ii], inputScalars.Nz[ii]);
#endif // END CUDA
				if (DEBUG) {
					mexPrintBase("global[0] = %u\n", global[0]);
					mexPrintBase("localBP[0] = %u\n", localBP[0]);
					mexPrintBase("localBP[1] = %u\n", localBP[1]);
					mexPrintBase("global[1] = %u\n", global[1]);
					mexPrintBase("global[2] = %u\n", global[2]);
					mexPrintBase("erotusBP[0] = %u\n", erotusBP[0][ii]);
					mexPrintBase("erotusBP[1] = %u\n", erotusBP[1][ii]);
#if defined(CUDA) || defined(HIP)
					mexPrintBase("kernelIndBPSubIter = %u\n", BPArgs.size());
#elif defined(OPENCL)
					mexPrintBase("kernelIndBPSubIter = %u\n", kernelIndBPSubIter);
					mexPrintBase("global.dimensions() = %u\n", global.dimensions());
					mexPrintBase("localBP.dimensions() = %u\n", localBP.dimensions());
#endif // END CUDA
					mexPrintBase("m_size = %u\n", m_size);
					mexPrintBase("nRowsD = %u\n", inputScalars.nRowsD);
					mexPrintBase("nColsD = %u\n", inputScalars.nColsD);
					mexPrintBase("length[indD] = %u\n", length[indD]);
					mexPrintBase("listmode = %u\n", inputScalars.listmode);
					mexPrintBase("im_dim = %u\n", inputScalars.im_dim[ii]);
					mexPrintBase("no_norm = %u\n", no_norm);
					mexPrintBase("memSize = %u\n", memSize);
					if (inputScalars.BPType == 4) {
						mexPrintBase("dL = %f\n", inputScalars.dL);
						mexPrintBase("kerroin4[ii] = %f\n", w_vec.kerroin4[ii]);
					}
					mexEval();
				}
				if (inputScalars.offset && inputScalars.BPType != 7) {
					KARG_METAL_SLOT(kernelIndBPSubIter, 1);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_T[osa_iter]);
				}
				if (inputScalars.BPType == 5 || inputScalars.BPType == 4 || inputScalars.BPType == 7) {
					KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, d_N[ii]);
					KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, b[ii]);
					KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, d[ii]);
					if (inputScalars.BPType == 5) {
						KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, inputScalars.d_Scale[ii]);
						KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, inputScalars.dSizeBP);
					}
					else if (inputScalars.BPType == 4) {
						KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, w_vec.kerroin4[ii]);
					}
				}
				if (inputScalars.BPType == 4) {
					KARG_METAL_SLOT(kernelIndBPSubIter, 2);
					if (inputScalars.useBuffers) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_output);
					}
					else
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_inputImage);
					if (inputScalars.CT && inputScalars.DSC > 0.f) {
						KARG_METAL_SLOT(kernelIndBPSubIter, 3);
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_angle);
						KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, inputScalars.DSC);
					}
					KARG_METAL_SLOT(kernelIndBPSubIter, 4);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, vec_opencl.d_rhs_os[uu]);
					KARG_METAL_SLOT(kernelIndBPSubIter, 5);
					if (compSens) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_xFull[0]);
					}
					else
						if (!inputScalars.loadTOF && inputScalars.listmode > 0) {
							KARG(kTemp, kernelBP, kernelIndBPSubIter, d_x[0][0]);
						}
						else
							KARG(kTemp, kernelBP, kernelIndBPSubIter, d_x[timestep][osa_iter]);
					KARG_METAL_SLOT(kernelIndBPSubIter, 6);
					if (compSens) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_zFull[0]);
					}
					else
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_z[timestep][osa_iter]);
					KARG_METAL_SLOT(kernelIndBPSubIter, 7);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_Summ[ee]);
				}
				else {
					KARG_METAL_SLOT(kernelIndBPSubIter, 5);
					if (compSens) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_xFull[0]);
					}
					else
						if (!inputScalars.loadTOF && inputScalars.listmode > 0) {
							KARG(kTemp, kernelBP, kernelIndBPSubIter, d_x[0][0]);
						}
							else
								KARG(kTemp, kernelBP, kernelIndBPSubIter, d_x[timestep][osa_iter]);
					KARG_METAL_SLOT(kernelIndBPSubIter, 6);
					if (compSens) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_zFull[0]);
					}
					else
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_z[timestep][osa_iter]);
					// Only when the kernel was built with -DGEOM5
#if !defined(METAL)
					if (inputScalars.BPType == 5 && inputScalars.CT && inputScalars.listmode == 0)
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_geomProj5[osa_iter]);
#endif
					KARG_METAL_SLOT(kernelIndBPSubIter, 2);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_inputImage);
					KARG_METAL_SLOT(kernelIndBPSubIter, 4);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, vec_opencl.d_rhs_os[uu]);
					KARG_METAL_SLOT(kernelIndBPSubIter, 7);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_Summ[ee]);
					if (inputScalars.meanBP) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_meanBP);
					}
				}
				if (inputScalars.normalization_correction) {
					KARG_METAL_SLOT(kernelIndBPSubIter, 8);
					// TODO: Listmode normalization
					//if (inputScalars.listmode > 0 && inputScalars.indexBased)
					//	status = kernelBP.setArg(kernelIndBPSubIter++, d_norm[0]);
					//else
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_norm[timestep][osa_iter]);
				}
			}
			else {
				if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0) {
					SET_LAUNCH_RANGE3(global,
						inputScalars.nRowsD + erotus[0],
						inputScalars.nColsD + erotus[1],
						length[indD],
						localBP);
				}
				else if (inputScalars.listmode > 0 && compSens) {
					SET_LAUNCH_RANGE3(global, static_cast<size_t>(inputScalars.det_per_ring) + erotusSens[0],
						static_cast<size_t>(inputScalars.det_per_ring) + erotusSens[1], 
						static_cast<size_t>(inputScalars.rings) * static_cast<size_t>(inputScalars.rings), localBP);
				}
				else {
					erotus[0] = length[indD] % local_size[0];
					if (erotus[0] > 0)
						erotus[0] = (local_size[0] - erotus[0]);
					// See the matching comment in forwardProjection: this launch is measurement driven and its
					// global range is 1D, so the 2D member range (forced to {16, 16} by fastPDHG for
					// projector types 4 and 5) has to be flattened, or global[1] == 1 is indivisible by it.
					SET_RANGE2(localBP, local_size[0], 1);
					SET_LAUNCH_RANGE3(global, length[indD] + erotus[0], 1, 1, localBP);
				}
				if (DEBUG) {
					mexPrintBase("global[0] = %u\n", global[0]);
					mexPrintBase("localBP[0] = %u\n", localBP[0]);
					mexPrintBase("localBP[1] = %u\n", localBP[1]);
					mexPrintBase("global[1] = %u\n", global[1]);
					mexPrintBase("global[2] = %u\n", global[2]);
					if (compSens) {
						mexPrintBase("erotusSens[0] = %u\n", erotusSens[0]);
						mexPrintBase("erotusSens[1] = %u\n", erotusSens[1]);
					}
					else {
						mexPrintBase("erotus[0] = %u\n", erotus[0]);
						mexPrintBase("erotus[1] = %u\n", erotus[1]);
					}
#if defined(CUDA) || defined(HIP)
					mexPrintBase("kernelIndBPSubIter = %u\n", BPArgs.size());
#elif defined(OPENCL)
					mexPrintBase("global.dimensions() = %u\n", global.dimensions());
					mexPrintBase("localBP.dimensions() = %u\n", localBP.dimensions());
					mexPrintBase("kernelIndBPSubIter = %u\n", kernelIndBPSubIter);
#endif // END CUDA
					mexPrintBase("m_size = %u\n", m_size);
					mexPrintBase("nRowsD = %u\n", inputScalars.nRowsD);
					mexPrintBase("nColsD = %u\n", inputScalars.nColsD);
					mexPrintBase("length[indD] = %u\n", length[indD]);
					mexPrintBase("listmode = %u\n", inputScalars.listmode);
					mexPrintBase("im_dim = %u\n", inputScalars.im_dim[ii]);
					mexPrintBase("no_norm = %u\n", no_norm);
					mexPrintBase("osa_iter = %u\n", osa_iter);
					mexPrintBase("memSize = %u\n", memSize);
					mexEval();
				}
				KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, d_N[ii]);
				KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, b[ii]);
				KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, bmax[ii]);
				KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, inputScalars.d_Scale4[ii]);
				if (inputScalars.largeDim) {
					// Non-CT BP 4 needs the same volume information as FP
					bzGlobalBP[0] = (ii == 0) ? inputScalars.lDimStruct.bz[0] : VEC_Z(b[ii]);
					bzGlobalBP[1] = (ii == 0) ? inputScalars.lDimStruct.bmaxZ[inputScalars.subsets - 1] : VEC_Z(bmax[ii]);
					KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, bzGlobalBP[0]);
					KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, bzGlobalBP[1]);
				}
				KARG_METAL_SLOT(kernelIndBPSubIter, 1);
				KARG(kTemp, kernelBP, kernelIndBPSubIter, d_output);
				KARG(kTemp, kernelBP, kernelIndBPSubIter, vec_opencl.d_rhs_os[uu]);
				KARG_METAL_SLOT(kernelIndBPSubIter, 3);
				if (compSens) {
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_xFull[0]);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_zFull[0]);
					KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, inputScalars.rings);
					KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, inputScalars.det_per_ring);
				}
				else {
					if (((inputScalars.listmode == 0 || inputScalars.indexBased) && !(inputScalars.CT || inputScalars.SPECT)) || (!inputScalars.loadTOF && inputScalars.listmode > 0)) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_x[0][0]);
					}
					else
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_x[timestep][osa_iter]);
					if ((inputScalars.CT || inputScalars.PET || inputScalars.SPECT)) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_z[timestep][osa_iter]);
					}
					else if (inputScalars.listmode > 0) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_z[0][0]);
					}
					else
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_z[timestep][inputScalars.osa_iter0]);
				}
				if (inputScalars.maskFP || inputScalars.maskBP) {
					if (inputScalars.maskFP) {
						KARG_METAL_SLOT(kernelIndBPSubIter, 7);
						if (inputScalars.useBuffers) {
							int subset = 0;
							if (inputScalars.maskFPZ > 1 && !(inputScalars.SPECT && inputScalars.maskFPZ == inputScalars.nHeads))
								subset = osa_iter;
							KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskFPB[timestep][subset]);
						}
						else
							if (inputScalars.maskFPZ > 1) {
								KARG(kTemp, kernelBP, kernelIndBPSubIter, (inputScalars.SPECT && inputScalars.maskFPZ == inputScalars.nHeads) ? d_maskFP3[timestep][0] : d_maskFP3[timestep][osa_iter]);
							}
							else
								KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskFP);
					}
					if (inputScalars.maskBP) {
						KARG_METAL_SLOT(kernelIndBPSubIter, 8);
						if (inputScalars.useBuffers) {
							KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskBPB[timestep][ii]);
						}
						else {
							if (inputScalars.maskBPZ > 1) {
								KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskBP3[timestep][ii]);
							}
							else
								KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskBP[timestep][ii]);
						}
					}
				}
				KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, length[indD]);
				if ((inputScalars.subsetType == 3 || inputScalars.subsetType == 6 || inputScalars.subsetType == 7) && inputScalars.subsets > 1 && !inputScalars.CT && !inputScalars.SPECT && !inputScalars.PET && inputScalars.listmode == 0) {
					KARG_METAL_SLOT(kernelIndBPSubIter, 9);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_xyindex[osa_iter]);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_zindex[osa_iter]);
				}
				if (inputScalars.listmode > 0 && inputScalars.indexBased && !compSens) {
					KARG_METAL_SLOT(kernelIndBPSubIter, 9);
					if (!inputScalars.loadTOF) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_trIndex[0][0]);
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_axIndex[0][0]);
					}
					else {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_trIndex[timestep][osa_iter]);
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_axIndex[timestep][osa_iter]);
					}
				}
				if (inputScalars.listmode > 0 && inputScalars.TOF) {
					KARG_METAL_SLOT(kernelIndBPSubIter, 11);
					if (!inputScalars.loadTOF) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_TOFIndex[0][0]);
					}
					else {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_TOFIndex[timestep][osa_iter]);
					}
				}
				if (inputScalars.raw) {
					KARG_METAL_SLOT(kernelIndBPSubIter, 12);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_L[osa_iter]);
					KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, inputScalars.det_per_ring);
				}
				if (inputScalars.normalization_correction) {
					KARG_METAL_SLOT(kernelIndBPSubIter, 13);
					//if (inputScalars.listmode > 0 && inputScalars.indexBased)
					//	status = kernelBP.setArg(kernelIndBPSubIter++, d_norm[0]);
					//else
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_norm[timestep][osa_iter]);
				}
				if (inputScalars.scatter) {
					KARG_METAL_SLOT(kernelIndBPSubIter, 14);
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_scat[timestep][osa_iter]);
				}
				KARG_METAL_SLOT(kernelIndBPSubIter, 15);
				KARG(kTemp, kernelBP, kernelIndBPSubIter, d_Summ[ee]);
				}
			KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, no_norm);
			if (inputScalars.CT && inputScalars.maskBP && (inputScalars.BPType == 4 || inputScalars.BPType == 5 || inputScalars.BPType == 7)) {
				KARG_METAL_SLOT(kernelIndBPSubIter, 9);
				if (inputScalars.useBuffers) {
					KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskBPB[timestep][ii]);
				}
				else {
					if (inputScalars.maskBPZ > 1) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskBP3[timestep][ii]);
					}
					else
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_maskBP[timestep][ii]);
				}
			}
			if (inputScalars.CT) {
				KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, length[indD]);
			}
			else {
				KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, m_size);
				KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, osa_iter);
			}
			KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, ii);
#ifndef METAL
			// fastPDHG backprojection inputs
			if (inputScalars.fastPDHG && inputScalars.CT && (inputScalars.BPType == 4 || inputScalars.BPType == 5)) {
				// Input the PDHG backprojection update only in the case of PDHG type algorithm
				// Add preconditioner if used
				const bool fastPrecondActive = !MethodList.BSREM && (w_vec.precondTypeIm[0] || w_vec.precondTypeIm[1]);
				if (fastStep != 0) {
					if (fastAlg == 0)
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_U);
					if (fastPrecondActive)
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_precond);
				}
				else {
					if (fastAlg == 0)
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_output);
					if (fastPrecondActive)
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_output);
				}
				if (fastNLMUsed) {
					// Use either the same cached image texture as with FPTypes 1-4
					if (fastStep == 0) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_inputImage);
					}
					// Or use the copy for FPType 5 or largeDim
					else if (inputScalars.FPType == 5 || inputScalars.largeDim) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_inputI);
					}
					else {
#if defined(CUDA) || defined(HIP)
						KARG(kTemp, kernelBP, kernelIndBPSubIter, FPTexCachePrior[0]);
#else
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_imageCache[0]);
#endif // END CUDA
					}
					if (MethodList.NLM) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_gaussianNLM);
						if (w_vec.NLM_anatomical)
							KARG(kTemp, kernelBP, kernelIndBPSubIter, d_urefIm);
						for (int nn = 0; nn < 6; nn++)
							KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, fastNLM[nn]);
					}
					else if (MethodList.RDP) {
						if (w_vec.RDPLargeNeighbor)
							KARG(kTemp, kernelBP, kernelIndBPSubIter, d_weights);
						KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, w_vec.RDP_gamma);
						if (w_vec.RDP_anatomical)
							// Anatomical reference image for RDP
							KARG(kTemp, kernelBP, kernelIndBPSubIter, d_RDPrefI);
					}
					else if (MethodList.GGMRF) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_weights);
						KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, w_vec.GGMRF_p);
						KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, w_vec.GGMRF_q);
						KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, w_vec.GGMRF_c);
						KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, w_vec.GGMRF_pqc);
					}
					else if (MethodList.hyperbolic) {
						KARG(kTemp, kernelBP, kernelIndBPSubIter, d_weights);
						KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, w_vec.data.SATVPhi);
					}
					else if (MethodList.TV) {
						KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, w_vec.data.SATVPhi);
					}
				}
				if (fastAlg == 0) {
					KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, fastTheta);
					KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, fastTau);
				}
				else {
					KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, fastLambda);
					KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, fastAlpha);
				}
				KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, fastBeta);
				KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, fastEpps);
				KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, fastPositivity);
				// largeDim offsets
				KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, fastImOffset);
				KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, fastPriorZOffset);
				KARG_SCALAR(kTemp, kernelBP, kernelIndBPSubIter, fastStep);
			}
#endif // END METAL
			}
		if (DEBUG || inputScalars.verbose >= 3)
			START_TIMER(tStart);
#if defined(CUDA) || defined(HIP)
		if (queueIdx > 0 && static_cast<size_t>(queueIdx) <= sideQueues.size()) {
			status = cuEventRecord(evMain, CLCommandQueue[0]);
			CUDA_CHECK(status, "Failed to record main stream event\n", -1);
			status = cuStreamWaitEvent(sideQueues[queueIdx - 1], evMain, 0);
			CUDA_CHECK(status, "Failed to make side stream wait for the main stream\n", -1);
			status = cuLaunchKernel(kernelBP, global[0], global[1], global[2], localBP[0], localBP[1], localBP[2], 0, sideQueues[queueIdx - 1], kTemp.data(), 0);
		}
		else
			status = cuLaunchKernel(kernelBP, global[0], global[1], global[2], localBP[0], localBP[1], localBP[2], 0, CLCommandQueue[0], kTemp.data(), 0);
		CUDA_CHECK(status, "Failed to launch backprojection kernel\n", -1);
#elif defined(METAL)
		{
			const MTL::Size threadsPerThreadgroup = MTL::Size::Make(localBP[0], localBP[1], localBP[2]);
			const MTL::Size threadgroupsPerGrid = MTL::Size::Make(global[0] / localBP[0], global[1] / localBP[1], global[2] / localBP[2]);
			encoder->dispatchThreadgroups(threadgroupsPerGrid, threadsPerThreadgroup);
			encoder->endEncoding();
			commandBuffer->commit();
			if (useMetalSideQueue) {
				sideCommandBuffers[queueIdx - 1] = commandBuffer;
			}
			else {
				commandBuffer->waitUntilCompleted();
				if (commandBuffer->status() == MTL::CommandBufferStatusError) {
					NS::Error* error = commandBuffer->error();
					const char* message = error && error->localizedDescription()
						? error->localizedDescription()->utf8String() : "unknown Metal error";
					mexPrintBase("Metal backprojection failed: %s\n", message);
					return -1;
				}
			}
		}
#elif defined(OPENCL)
		if (queueIdx > 0 && static_cast<size_t>(queueIdx) <= sideQueues.size()) {
			cl::Event evMain;
			status = CLCommandQueue[0].enqueueMarkerWithWaitList(nullptr, &evMain);
			OCL_CHECK(status, "Failed to enqueue main queue marker\n", -1);
			std::vector<cl::Event> waitList = { evMain };
			status = sideQueues[queueIdx - 1].enqueueBarrierWithWaitList(&waitList);
			OCL_CHECK(status, "Failed to enqueue side queue barrier\n", -1);
			status = sideQueues[queueIdx - 1].enqueueNDRangeKernel(kernelBP, cl::NDRange(), global, localBP, NULL);
			OCL_CHECK(status, "\n", -1);
			status = sideQueues[queueIdx - 1].flush();
			OCL_CHECK(status, "Failed to flush side queue\n", -1);
		}
		else {
			status = CLCommandQueue[0].enqueueNDRangeKernel(kernelBP, cl::NDRange(), global, localBP, NULL);
			OCL_CHECK(status, "\n", -1);
		}
#endif // END CUDA
		if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Backprojection kernel launched successfully\n");
		}
#if defined(CUDA) || defined(HIP)
		if (DEBUG || inputScalars.verbose >= 3)
			cuEventRecord(tEnd, queueIdx > 0 && static_cast<size_t>(queueIdx) <= sideQueues.size() ? sideQueues[queueIdx - 1] : CLCommandQueue[0]);
#ifndef AF
		status = cuCtxSynchronize();
		CUDA_CHECK(status, "Synchronization failed after backprojection\n", -1);
#endif
#elif defined(OPENCL)
#ifndef AF
		status = CLCommandQueue[0].finish();
		OCL_CHECK(status, "\n", -1);
#endif
#endif // END CUDA
		if (inputScalars.listmode > 0 && compSens) {
			kernelBP = kernelApu;
		}
		if (DEBUG || inputScalars.verbose >= 3) {
#if defined(CUDA) || defined(HIP)
			PRINT_TIMER(tStart, tEnd, "Backprojection completed in %f seconds\n");
			cuEventDestroy(tStart);
			cuEventDestroy(tEnd);
#elif defined(OPENCL)
			if (queueIdx > 0 && static_cast<size_t>(queueIdx) <= sideQueues.size())
				sideQueues[queueIdx - 1].finish();
			else
				CLCommandQueue[0].finish();
			STOP_TIMER(tEnd);
			PRINT_TIMER(tStart, tEnd, "Backprojection completed in %f seconds\n");
#elif defined(METAL)
			STOP_TIMER(tEnd);
			if (useMetalSideQueue)
				PRINT_TIMER(tStart, tEnd, "Backprojection submitted in %f seconds\n");
			else
				PRINT_TIMER(tStart, tEnd, "Backprojection completed in %f seconds\n");
#endif // END CUDA
		}
		return 0;
		}


	/// <summary>
	/// Release buffers needed only by the initial computation of the sensitivity image including all measurements (e.g. image-based preconditioners 2-3)
	/// </summary>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <returns></returns>
	inline void releaseBuffer(const scalarStruct & inputScalars) {
		d_xFull.clear();
		d_zFull.clear();
		if (inputScalars.size_norm > 1 && inputScalars.normalization_correction) {
			d_normFull.clear();
		}
		if (inputScalars.size_scat > 1 && inputScalars.scatter == 1U) {
			d_scatFull.clear();
		}
	}

	/// <summary>
	/// Get the total global memory of the selected device
	/// </summary>
	inline int64_t getGlobalMem() {
		STATUS_t status = SUCCESS_VALUE;
		uint64_t mem = 0ULL;
		uint64_t memFree = 0ULL;
		uint64_t memLoc = 0ULL;
#if defined(CUDA) || defined(HIP)
		size_t memTotalCuda = 0ULL;
		size_t memFreeCuda = 0ULL;
		status = cuMemGetInfo(&memFreeCuda, &memTotalCuda);
		CUDA_CHECK(status, "\n", -1);
		mem = static_cast<uint64_t>(memTotalCuda);
		memFree = static_cast<uint64_t>(memFreeCuda);
		int sharedMem = 0;
		cuDeviceGetAttribute(&sharedMem, CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK, CUDeviceID[0]);
		memLoc = static_cast<uint64_t>(sharedMem);
#elif defined(METAL)
		if (mtlDevice)
			mem = static_cast<uint64_t>(mtlDevice->recommendedMaxWorkingSetSize());
		else
			status = -1;
#elif defined(OPENCL)
		mem = static_cast<uint64_t>(CLDeviceID[0].getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>(&status));
		OCL_CHECK(status, "\n", -1);
		memLoc = static_cast<uint64_t>(CLDeviceID[0].getInfo<CL_DEVICE_LOCAL_MEM_SIZE>(&status));
#endif // END CUDA
		CHECK(status, "\n", -1);
		if (DEBUG) {
			mexPrintBase("mem_loc = %u\n", memLoc);
			if (memFree > 0ULL)
				mexPrintBase("memFree = %u\n", memFree);
		}
		return static_cast<int64_t>(mem);
	}

#if !defined(METAL)
	/// <summary>
	/// Compute median root prior (MRP)
	/// </summary>
	/// <param name="padd the padded input array (current estimate)"></param>
	/// <param name="grad the output gradient array"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <returns></returns>
#if defined(CUDA) || defined(HIP)
	inline int computeMRP(const scalarStruct & inputScalars, const uint64_t global_size[]) {
		std::vector<void*> kArgs;
#elif defined(OPENCL)
	inline int computeMRP(const scalarStruct & inputScalars, const uint64_t gSize[]) {
		UINT32_t kernelIndMed = 0U;
		uint64_t erotus[2] = { gSize[0] % localPrior[0], gSize[1] % localPrior[1] };
		cl::NDRange global_size(gSize[0] + (localPrior[0] - erotus[0]), gSize[1] + (localPrior[1] - erotus[1]), gSize[2]);
#endif // END CUDA
		TimerPoint tStart, tEnd;
		STATUS_t status = SUCCESS_VALUE;
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrint("Starting " BACKEND_STR " median kernel computation");
		if (DEBUG || inputScalars.verbose >= 3) {
			INIT_TIMER(tStart, tEnd);
		}
#if defined(CUDA) || defined(HIP)
		unsigned int gSize[3];
		unsigned int erotus[2];
		erotus[0] = localPrior[0] - (global_size[0] % localPrior[0]);
		erotus[1] = localPrior[1] - (global_size[1] % localPrior[1]);
		gSize[0] = (global_size[0] + erotus[0]) / localPrior[0];
		gSize[1] = (global_size[1] + erotus[1]) / localPrior[1];
		gSize[2] = global_size[2];
#endif // END CUDA
		//FINISH_QUEUE(status, "\n", -1);
		if (DEBUG) {
			mexPrintBase("global_size[0] = %d\n", global_size[0]);
			mexPrintBase("global_size[1] = %d\n", global_size[1]);
			mexPrintBase("global_size[2] = %d\n", global_size[2]);
			mexPrintBase("erotus[0] = %d\n", erotus[0]);
			mexPrintBase("erotus[1] = %d\n", erotus[1]);
			mexPrintBase("gSize[0] = %d\n", gSize[0]);
			mexPrintBase("gSize[1] = %d\n", gSize[1]);
			mexEval();
		}
		KARG(kArgs, kernelMed, kernelIndMed, d_inputB);
		KARG(kArgs, kernelMed, kernelIndMed, d_W);
		KARG(kArgs, kernelMed, kernelIndMed, d_N[0]);
		KARG(kArgs, kernelMed, kernelIndMed, d_NOrig);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.useBuffers) {
				KARG(kArgs, kernelMed, kernelIndMed, d_maskPriorB);
			}
			else
#if defined(OPENCL)
				if (inputScalars.maskBPZ > 1) {
					KARG(kArgs, kernelMed, kernelIndMed, d_maskPrior3);
				}
				else
#endif // END CUDA
					KARG(kArgs, kernelMed, kernelIndMed, d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			KARG(kArgs, kernelMed, kernelIndMed, d_eFOVIndices);
		if (DEBUG || inputScalars.verbose >= 3)
			START_TIMER(tStart);
#if defined(CUDA) || defined(HIP)
		status = cuLaunchKernel(kernelMed, gSize[0], gSize[1], gSize[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
#elif defined(OPENCL)
		status = CLCommandQueue[0].enqueueNDRangeKernel(kernelMed, cl::NullRange, global_size, localPrior);
#endif // END CUDA
		CHECK(status, "Failed to launch the Median filter kernel\n", -1);
		if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Median kernel launched successfully\n");
		}
		if (DEBUG || inputScalars.verbose >= 3) {
			FINISH_QUEUE(status, "Queue finish failed after MRP kernel\n", status);
			STOP_TIMER(tEnd);
			PRINT_TIMER(tStart, tEnd, BACKEND_STR " MRP kernel completed in %f seconds\n");
		}
		return 0;
	}

	/// <summary>
	/// Non-local means (NLM) prior
	/// </summary>
	/// <param name="grad the output gradient array"></param>
	/// <param name="im the input array (current estimate)"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="w_vec specifies some of the special options/parameters used"></param>
	/// <returns></returns>
#if defined(CUDA) || defined(HIP)
	inline int computeNLM(const scalarStruct & inputScalars, Weighting & w_vec, float beta, const int kk = 0) {
#elif defined(OPENCL)
	inline int computeNLM(const scalarStruct & inputScalars, Weighting & w_vec, const float beta, const int kk = 0) {
#endif // END CUDA
		TimerPoint tStart, tEnd;
		if (DEBUG || inputScalars.verbose >= 3) {
			INIT_TIMER(tStart, tEnd);
		}
		STATUS_t status = SUCCESS_VALUE;
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrint("Starting " BACKEND_STR " NLM gradient computation");
#if defined(CUDA) || defined(HIP)
		std::vector<void*> kArgs;
		float apu = inputScalars.epps;
#elif defined(OPENCL)
		const INT3_t searchWindow = { static_cast<INT32_t>(w_vec.Ndx) , static_cast<INT32_t>(w_vec.Ndy) , static_cast<INT32_t>(w_vec.Ndz) };
		const INT3_t patchWindow = { static_cast<INT32_t>(w_vec.Nlx) , static_cast<INT32_t>(w_vec.Nly) , static_cast<INT32_t>(w_vec.Nlz) };
		UINT32_t kernelIndNLM = 0U;
#endif // END CUDA
		uint32_t Nz, NzOrig;
		UINT2_t nOffset;
		if (inputScalars.largeDim)
			Nz = inputScalars.lDimStruct.NzPr[kk];
		else
			Nz = inputScalars.Nz[0];
		if (inputScalars.largeDim) {
			SET_RANGE_Z(globalPrior, Nz);
			NzOrig = VEC_Z(d_N[0]);
			VEC_Z(d_N[0]) = Nz;
		}
		if (kk == 0 && inputScalars.largeDim)
			nOffset = { 0, Nz - w_vec.Nlz - w_vec.Ndz };
		else if (kk < inputScalars.subsetsUsed - 1 && kk > 0 && inputScalars.largeDim)
			nOffset = { w_vec.Nlz + w_vec.Ndz, Nz - w_vec.Nlz - w_vec.Ndz };
		else if (inputScalars.largeDim)
			nOffset = { w_vec.Nlz + w_vec.Ndz, Nz };
		if (DEBUG) {
			mexPrintBase("w_vec.Ndx = %u\n", w_vec.Ndx);
			mexPrintBase("w_vec.Ndy = %u\n", w_vec.Ndy);
			mexPrintBase("w_vec.Ndz = %u\n", w_vec.Ndz);
			mexPrintBase("w_vec.Nlx = %u\n", w_vec.Nlx);
			mexPrintBase("w_vec.Nly = %u\n", w_vec.Nly);
			mexPrintBase("w_vec.Nlz = %u\n", w_vec.Nlz);
			mexPrintBase("erotusPrior[0] = %u\n", erotusPrior[0]);
			mexPrintBase("erotusPrior[1] = %u\n", erotusPrior[1]);
			mexPrintBase("erotusPrior[2] = %u\n", erotusPrior[2]);
			mexPrintBase("globalPrior[0] = %u\n", globalPrior[0]);
			mexPrintBase("globalPrior[1] = %u\n", globalPrior[1]);
			mexPrintBase("globalPrior[2] = %u\n", globalPrior[2]);
			mexPrintBase("localPrior[0] = %u\n", localPrior[0]);
			mexPrintBase("localPrior[1] = %u\n", localPrior[1]);
			mexPrintBase("localPrior[2] = %u\n", localPrior[2]);
			mexPrintBase("w_vec.h2 = %f\n", w_vec.h2);
			mexPrintBase("Nz = %u\n", Nz);
			mexPrintBase("kk = %u\n", kk);
			mexPrintBase("nOffset.x = %u\n", VEC_X(nOffset));
			mexPrintBase("nOffset.y = %u\n", VEC_Y(nOffset));
			mexPrintBase("d_N[0].z = %u\n", VEC_Z(d_N[0]));
			mexPrintBase("w_vec.RDP_gamma = %f\n", w_vec.RDP_gamma);
			mexPrintBase("useImages = %d\n", inputScalars.useImages);
			mexEval();
		}
		KARG(kArgs, kernelNLM, kernelIndNLM, d_W);
		if (inputScalars.useImages) {
			// FPType == 5 caches integral images rather than the estimate, so the prior uses its own copy
			if (inputScalars.FPType == 5) {
				KARG(kArgs, kernelNLM, kernelIndNLM, d_inputI);
			}
			else {
#if defined(CUDA) || defined(HIP)
				KARG(kArgs, kernelNLM, kernelIndNLM, FPTexCachePrior[kk]);
#else
				KARG(kArgs, kernelNLM, kernelIndNLM, d_imageCache[kk]);
#endif
			}
		}
		else {
			KARG(kArgs, kernelNLM, kernelIndNLM, d_inputB);
		}
		KARG(kArgs, kernelNLM, kernelIndNLM, d_gaussianNLM);
		KARG(kArgs, kernelNLM, kernelIndNLM, d_N[0]);
		KARG(kArgs, kernelNLM, kernelIndNLM, d_NOrig);
		KARG(kArgs, kernelNLM, kernelIndNLM, w_vec.h2);
#if defined(CUDA) || defined(HIP)
		KARG(kArgs, kernelNLM, kernelIndNLM, apu);
#elif defined(OPENCL)
		KARG(kArgs, kernelNLM, kernelIndNLM, inputScalars.epps);
#endif // END CUDA
		KARG(kArgs, kernelNLM, kernelIndNLM, beta);
		if (w_vec.NLRD || w_vec.NLLange || w_vec.NLGGMRF)
			KARG(kArgs, kernelNLM, kernelIndNLM, w_vec.RDP_gamma);
		if (w_vec.NLGGMRF) {
			KARG(kArgs, kernelNLM, kernelIndNLM, w_vec.GGMRF_p);
			KARG(kArgs, kernelNLM, kernelIndNLM, w_vec.GGMRF_q);
			KARG(kArgs, kernelNLM, kernelIndNLM, w_vec.GGMRF_c);
		}
		if (w_vec.NLAdaptive)
			KARG(kArgs, kernelNLM, kernelIndNLM, w_vec.NLAdaptiveConstant);
		if (w_vec.NLM_anatomical)
			if (inputScalars.useImages) {
				KARG(kArgs, kernelNLM, kernelIndNLM, d_urefIm);
			}
			else {
				KARG(kArgs, kernelNLM, kernelIndNLM, d_uref);
			}
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.useBuffers) {
				KARG(kArgs, kernelNLM, kernelIndNLM, d_maskPriorB);
			}
			else
#if defined(OPENCL)
				if (inputScalars.maskBPZ > 1) {
					KARG(kArgs, kernelNLM, kernelIndNLM, d_maskPrior3);
				}
				else
#endif // END CUDA
					KARG(kArgs, kernelNLM, kernelIndNLM, d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			KARG(kArgs, kernelNLM, kernelIndNLM, d_eFOVIndices);
		if (inputScalars.largeDim)
			KARG(kArgs, kernelNLM, kernelIndNLM, nOffset);
		//Compute the kernel
		if (DEBUG || inputScalars.verbose >= 3)
			START_TIMER(tStart);
#if defined(CUDA) || defined(HIP)
		status = cuLaunchKernel(kernelNLM, globalPrior[0], globalPrior[1], globalPrior[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		CUDA_CHECK(status, "Failed to launch the NLM kernel\n", status);

		if (inputScalars.useImages && inputScalars.FPType == 5) {
			status = cuTexObjectDestroy(d_inputI);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
			}
			status = cuArrayDestroy(imArray);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
			}
		}
#elif defined(OPENCL)
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelNLM, cl::NullRange, globalPrior, localPrior);
		OCL_CHECK(status, "Failed to launch the NLM kernel\n", -1);
#endif // END CUDA
		if (DEBUG || inputScalars.verbose >= 3) {
			FINISH_QUEUE(status, "Queue finish failed after NLM kernel\n", status);
			STOP_TIMER(tEnd);
			PRINT_TIMER(tStart, tEnd, BACKEND_STR " NLM gradient completed in %f seconds\n");
		}
		if (inputScalars.largeDim)
			VEC_Z(d_N[0]) = NzOrig;
		return 0;
	}

	/// <summary>
	/// Compute relative difference prior (RDP)
	/// </summary>
	/// <param name="grad the output gradient array"></param>
	/// <param name="im the input array (current estimate)"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="gamma controls the shape of the prior"></param>
	/// <param name="weights_RDP (UNUSED) the voxel weights for RDP"></param>
	/// <returns></returns>
#if defined(CUDA) || defined(HIP)
	inline int computeRDP(const scalarStruct & inputScalars, float gamma, const Weighting & w_vec, float beta, const int kk = 0, const bool RDPLargeNeighbor = false, const bool useRDPRef = false) {
#elif defined(OPENCL)
	inline int computeRDP(const scalarStruct & inputScalars, const float gamma, const Weighting & w_vec, const float beta, const int kk = 0, const bool RDPLargeNeighbor = false, const bool useRDPRef = false) {
#endif // END CUDA
		STATUS_t status = SUCCESS_VALUE;
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrint("Starting " BACKEND_STR " RDP gradient computation");
#if defined(CUDA) || defined(HIP)
		std::vector<void*> kArgs;
		float apu = inputScalars.epps;
#endif // END CUDA
		TimerPoint tStart, tEnd;
		if (DEBUG || inputScalars.verbose >= 3) {
			INIT_TIMER(tStart, tEnd);
		}
#if defined(OPENCL)
		UINT32_t kernelIndRDP = 0U;
		if (inputScalars.largeDim)
			SET_RANGE_Z(globalPrior, inputScalars.Nz[0]);
#endif // END CUDA
		uint32_t Nz, NzOrig;
		UINT2_t nOffset;
		if (inputScalars.largeDim)
			Nz = inputScalars.lDimStruct.NzPr[kk];
		else
			Nz = inputScalars.Nz[0];
		if (inputScalars.largeDim) {
			SET_RANGE_Z(globalPrior, Nz);
			NzOrig = VEC_Z(d_N[0]);
			VEC_Z(d_N[0]) = Nz;
		}
		if (kk == 0 && inputScalars.largeDim)
			nOffset = { 0, NzOrig };
		else if (kk < inputScalars.subsetsUsed - 1 && kk > 0 && inputScalars.largeDim)
			nOffset = { (Nz - NzOrig) / 2, (Nz + NzOrig) / 2 };
		else if (inputScalars.largeDim)
			nOffset = { Nz - NzOrig, Nz };
		//FINISH_QUEUE(status, "Queue finish failed before RDP kernel\n", -1);
		if (DEBUG) {
			mexPrintBase("inputScalars.epps = %.9f\n", inputScalars.epps);
			mexPrintBase("gamma = %f\n", gamma);
			mexPrintBase("inputScalars.Nx = %d\n", inputScalars.Nx[0]);
			mexPrintBase("inputScalars.Ny = %d\n", inputScalars.Ny[0]);
			mexPrintBase("inputScalars.Nz * inputScalars.nRekos = %d\n", inputScalars.Nz[0] * inputScalars.nRekos);
			mexPrintBase("globalPrior[0] = %d\n", globalPrior[0]);
			mexPrintBase("globalPrior[1] = %d\n", globalPrior[1]);
			mexPrintBase("globalPrior[2] = %d\n", globalPrior[2]);
			mexPrintBase("RDPLargeNeighbor = %d\n", RDPLargeNeighbor);
			mexEval();
		}
		KARG(kArgs, kernelRDP, kernelIndRDP, d_W);
		if (inputScalars.useImages) {
			// FPType == 5 caches integral images rather than the estimate, so the prior uses its own copy
			if (inputScalars.FPType == 5) {
				KARG(kArgs, kernelRDP, kernelIndRDP, d_inputI);
			}
			else {
#if defined(CUDA) || defined(HIP)
				KARG(kArgs, kernelRDP, kernelIndRDP, FPTexCachePrior[kk]);
#else
				KARG(kArgs, kernelRDP, kernelIndRDP, d_imageCache[kk]);
#endif
			}
		}
		else {
			KARG(kArgs, kernelRDP, kernelIndRDP, d_inputB);
		}
		KARG(kArgs, kernelRDP, kernelIndRDP, d_N[0]);
		KARG(kArgs, kernelRDP, kernelIndRDP, d_NOrig);
		KARG(kArgs, kernelRDP, kernelIndRDP, gamma);
#if defined(CUDA) || defined(HIP)
		KARG(kArgs, kernelRDP, kernelIndRDP, apu);
#elif defined(OPENCL)
		KARG(kArgs, kernelRDP, kernelIndRDP, inputScalars.epps);
#endif // END CUDA
		KARG(kArgs, kernelRDP, kernelIndRDP, beta);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.useBuffers) {
				KARG(kArgs, kernelRDP, kernelIndRDP, d_maskPriorB);
			}
			else
#if defined(OPENCL)
				if (inputScalars.maskBPZ > 1) {
					KARG(kArgs, kernelRDP, kernelIndRDP, d_maskPrior3);
				}
				else
#endif // END CUDA
					KARG(kArgs, kernelRDP, kernelIndRDP, d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			KARG(kArgs, kernelRDP, kernelIndRDP, d_eFOVIndices);
		if (RDPLargeNeighbor) {
			KARG(kArgs, kernelRDP, kernelIndRDP, d_weights);
			if (useRDPRef)
				if (inputScalars.useImages) {
					KARG(kArgs, kernelRDP, kernelIndRDP, d_RDPrefI);
				}
				else {
					KARG(kArgs, kernelRDP, kernelIndRDP, d_RDPref);
				}
		}
		if (inputScalars.largeDim)
			KARG(kArgs, kernelRDP, kernelIndRDP, nOffset);
		// Compute the kernel
		if (DEBUG || inputScalars.verbose >= 3)
			START_TIMER(tStart);
#if defined(CUDA) || defined(HIP)
		status = cuLaunchKernel(kernelRDP, globalPrior[0], globalPrior[1], globalPrior[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		CUDA_CHECK(status, "Failed to launch the RDP kernel\n", -1);

		if (inputScalars.useImages) {
			if (inputScalars.FPType == 5) {
				status = cuTexObjectDestroy(d_inputI);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
				}
				status = cuArrayDestroy(imArray);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
				}
			}
			if (RDPLargeNeighbor && useRDPRef) {
				status = cuTexObjectDestroy(d_RDPrefI);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
				}
			}
		}
#elif defined(OPENCL)
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelRDP, cl::NullRange, globalPrior, localPrior);
		OCL_CHECK(status, "Failed to launch the RDP kernel\n", -1);
#endif // END CUDA
		if (DEBUG || inputScalars.verbose >= 3) {
			FINISH_QUEUE(status, "Queue finish failed after RDP kernel\n", -1);
			STOP_TIMER(tEnd);
			PRINT_TIMER(tStart, tEnd, BACKEND_STR " RDP gradient completed in %f seconds\n");
		}
		if (inputScalars.largeDim)
			VEC_Z(d_N[0]) = NzOrig;
		return 0;
	}

	/// <summary>
	/// Compute relative generalized Gaussian Markov random field prior (GGMRF)
	/// </summary>
	/// <param name="grad the output gradient array"></param>
	/// <param name="im the input array (current estimate)"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="p constant controlling the powers near from the origin"></param>
	/// <param name="q constant controlling the powers distant from the origin"></param>
	/// <param name="c constant controlling the approximate threshold of transition between low and high contrast regions"></param>
	/// <param name="beta regularization parameter"></param>
	/// <returns></returns>
#if defined(CUDA) || defined(HIP)
	inline int computeGGMRF(const scalarStruct & inputScalars, float p, float q, float c, float pqc, const Weighting & w_vec, float beta, const int kk = 0) {
#elif defined(OPENCL)
	inline int computeGGMRF(const scalarStruct & inputScalars, const float p, const float q, const float c, const float pqc, const Weighting & w_vec, const float beta, const int kk = 0) {
#endif // END CUDA
		STATUS_t status = SUCCESS_VALUE;
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrint("Starting " BACKEND_STR " GGMRF gradient computation");
		TimerPoint tStart, tEnd;
		if (DEBUG || inputScalars.verbose >= 3) {
			INIT_TIMER(tStart, tEnd);
		}
#if defined(OPENCL)
		if (inputScalars.largeDim)
			SET_RANGE_Z(globalPrior, inputScalars.Nz[0]);
#endif // END CUDA
		uint32_t Nz, NzOrig;
		UINT2_t nOffset;
		if (inputScalars.largeDim)
			Nz = inputScalars.lDimStruct.NzPr[kk];
		else
			Nz = inputScalars.Nz[0];
		if (inputScalars.largeDim) {
			SET_RANGE_Z(globalPrior, Nz);
			NzOrig = VEC_Z(d_N[0]);
			VEC_Z(d_N[0]) = Nz;
		}
		if (kk == 0 && inputScalars.largeDim)
			nOffset = { 0, Nz - w_vec.Ndz };
		else if (kk < inputScalars.subsetsUsed - 1 && kk > 0 && inputScalars.largeDim)
			nOffset = { w_vec.Ndz, Nz - w_vec.Ndz };
		else if (inputScalars.largeDim)
			nOffset = { w_vec.Ndz, Nz };
#if defined(CUDA) || defined(HIP)
		std::vector<void*> kArgs;
#elif defined(OPENCL)
		UINT32_t kernelIndGGMRF = 0U;
#endif // END CUDA
		if (DEBUG) {
			mexPrintBase("p = %f\n", p);
			mexPrintBase("q = %f\n", q);
			mexPrintBase("c = %f\n", c);
			mexPrintBase("pqc = %f\n", pqc);
			mexPrintBase("inputScalars.Nx = %d\n", inputScalars.Nx[0]);
			mexPrintBase("inputScalars.Ny = %d\n", inputScalars.Ny[0]);
			mexPrintBase("inputScalars.Nz * inputScalars.nRekos = %d\n", inputScalars.Nz[0] * inputScalars.nRekos);
			mexPrintBase("globalPrior[0] = %d\n", globalPrior[0]);
			mexPrintBase("globalPrior[1] = %d\n", globalPrior[1]);
			mexPrintBase("globalPrior[2] = %d\n", globalPrior[2]);
			mexEval();
		}
		KARG(kArgs, kernelGGMRF, kernelIndGGMRF, d_W);
		if (inputScalars.useImages) {
			// FPType == 5 caches integral images rather than the estimate, so the prior uses its own copy
			if (inputScalars.FPType == 5) {
				KARG(kArgs, kernelGGMRF, kernelIndGGMRF, d_inputI);
			}
			else {
#if defined(CUDA) || defined(HIP)
				KARG(kArgs, kernelGGMRF, kernelIndGGMRF, FPTexCachePrior[kk]);
#else
				KARG(kArgs, kernelGGMRF, kernelIndGGMRF, d_imageCache[kk]);
#endif
			}
		}
		else {
			KARG(kArgs, kernelGGMRF, kernelIndGGMRF, d_inputB);
		}
		KARG(kArgs, kernelGGMRF, kernelIndGGMRF, d_weights);
		KARG(kArgs, kernelGGMRF, kernelIndGGMRF, d_N[0]);
		KARG(kArgs, kernelGGMRF, kernelIndGGMRF, p);
		KARG(kArgs, kernelGGMRF, kernelIndGGMRF, q);
		KARG(kArgs, kernelGGMRF, kernelIndGGMRF, c);
		KARG(kArgs, kernelGGMRF, kernelIndGGMRF, pqc);
		KARG(kArgs, kernelGGMRF, kernelIndGGMRF, beta);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.useBuffers) {
				KARG(kArgs, kernelGGMRF, kernelIndGGMRF, d_maskPriorB);
			}
			else
#if defined(OPENCL)
				if (inputScalars.maskBPZ > 1) {
					KARG(kArgs, kernelGGMRF, kernelIndGGMRF, d_maskPrior3);
				}
				else
#endif // END CUDA
					KARG(kArgs, kernelGGMRF, kernelIndGGMRF, d_maskPrior);
		if (inputScalars.largeDim)
			KARG(kArgs, kernelGGMRF, kernelIndGGMRF, nOffset);
		// Compute the kernel
		if (DEBUG || inputScalars.verbose >= 3)
			START_TIMER(tStart);
#if defined(CUDA) || defined(HIP)
		status = cuLaunchKernel(kernelGGMRF, globalPrior[0], globalPrior[1], globalPrior[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		CUDA_CHECK(status, "Failed to launch the GGMRF kernel\n", -1);

		if (inputScalars.useImages && inputScalars.FPType == 5) {
			status = cuTexObjectDestroy(d_inputI);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
			}
			status = cuArrayDestroy(imArray);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
			}
		}
#elif defined(OPENCL)
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelGGMRF, cl::NullRange, globalPrior, localPrior);
		OCL_CHECK(status, "Failed to launch the GGMRF kernel\n", -1);
#endif // END CUDA
		if (DEBUG || inputScalars.verbose >= 3) {
			FINISH_QUEUE(status, "Queue finish failed after GGMRF kernel\n", -1);
			STOP_TIMER(tEnd);
			PRINT_TIMER(tStart, tEnd, BACKEND_STR " GGMRF gradient completed in %f seconds\n");
		}
		if (inputScalars.largeDim)
			VEC_Z(d_N[0]) = NzOrig;
		return 0;
	}


#if defined(CUDA) || defined(HIP)
	inline int ProxHelperQ(float alpha, const uint64_t globalQ) {
		std::vector<void*> kArgs;
		UINT32_t kernelIndProxRDP = 0U;
#elif defined(OPENCL)
	inline int ProxHelperQ(const float alpha, const uint64_t gQ) {
		cl::NDRange globalQ = { static_cast<cl::size_type>(gQ) };
		UINT32_t kernelIndProxRDP = 0U;
#endif // END CUDA
		STATUS_t status = SUCCESS_VALUE;
		//FINISH_QUEUE(status, "Queue finish failed before proximal RDP helper kernel\n", -1);
		KARG(kArgs, kernelProxq, kernelIndProxRDP, d_qX);
		KARG(kArgs, kernelProxq, kernelIndProxRDP, alpha);
		// Compute the kernel
#if defined(CUDA) || defined(HIP)
		status = cuLaunchKernel(kernelProxq, globalQ / 64ULL, 1, 1, 64, 1, 1, 0, CLCommandQueue[0], kArgs.data(), NULL);
		CUDA_CHECK(status, "Failed to launch the Proximal RDP helper kernel\n", -1);
#elif defined(OPENCL)
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelProxq, cl::NullRange, globalQ, cl::NullRange);
		OCL_CHECK(status, "Failed to launch the Proximal RDP helper kernel\n", -1);
#endif // END CUDA
		//FINISH_QUEUE(status, "Queue finish failed after proximal RDP helper kernel\n", -1);
		return status;
	}

	/// <summary>
	/// Compute either the projection from an L1 or L2 ball for the TV prior
	/// </summary>
	/// <param name="q the input TV array"></param>
	/// <param name="alpha the regularization parameter"></param>
	/// <param name="L2Ball if true, computes the projection from an L2 ball, otherwise from the L1 ball"></param>
	/// <returns></returns>
#if defined(CUDA) || defined(HIP)
	inline int ProxTVHelperQ(float alpha, const uint64_t globalQ) {
		std::vector<void*> kArgs;
#elif defined(OPENCL)
	inline int ProxTVHelperQ(const float alpha, const uint64_t gQ) {
		cl::NDRange globalQ = { static_cast<cl::size_type>(gQ) };
		UINT32_t kernelIndCPTV = 0U;
#endif // END CUDA
		STATUS_t status = SUCCESS_VALUE;
		//FINISH_QUEUE(status, "Queue finish failed before proximal TV kernel\n", -1);
		KARG(kArgs, kernelProxTVq, kernelIndCPTV, d_qX);
		KARG(kArgs, kernelProxTVq, kernelIndCPTV, d_qY);
		KARG(kArgs, kernelProxTVq, kernelIndCPTV, d_qZ);
		KARG(kArgs, kernelProxTVq, kernelIndCPTV, alpha);
		// Compute the kernel
#if defined(CUDA) || defined(HIP)
		status = cuLaunchKernel(kernelProxTVq, globalQ / 64ULL, 1, 1, 64, 1, 1, 0, CLCommandQueue[0], kArgs.data(), NULL);
		CUDA_CHECK(status, "Failed to launch the Proximal TV kernel\n", -1);
#elif defined(OPENCL)
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelProxTVq, cl::NullRange, globalQ, cl::NullRange);
		OCL_CHECK(status, "Failed to launch the Proximal TV kernel\n", -1);
#endif // END CUDA
		//FINISH_QUEUE(status, "Queue finish failed after proximal TV kernel\n", -1);
		return status;
	}

	/// <summary>
	/// Compute either the projection from an L1 or L2 ball for the TGV prior
	/// </summary>
	/// <param name="q first half of the input TGV array"></param>
	/// <param name="q2 second half of the input TGV array"></param>
	/// <param name="alpha the regularization parameter"></param>
	/// <returns></returns>
#if defined(CUDA) || defined(HIP)
	inline int ProxTGVHelperQ(const scalarStruct & inputScalars, float alpha, const uint64_t globalQ) {
		std::vector<void*> kArgs;
#elif defined(OPENCL)
	inline int ProxTGVHelperQ(const scalarStruct & inputScalars, const float alpha, const uint64_t globalQ) {
		UINT32_t kernelIndCPTV = 0U;
#endif // END CUDA
		STATUS_t status = SUCCESS_VALUE;
		//FINISH_QUEUE(status, "Queue finish failed before proximal TGV kernel\n", -1);
		KARG(kArgs, kernelProxTGVq, kernelIndCPTV, d_rX);
		KARG(kArgs, kernelProxTGVq, kernelIndCPTV, d_rY);
		if (!inputScalars.TGV2D)
			KARG(kArgs, kernelProxTGVq, kernelIndCPTV, d_rZ);
		KARG(kArgs, kernelProxTGVq, kernelIndCPTV, d_rXY);
		if (!inputScalars.TGV2D) {
			KARG(kArgs, kernelProxTGVq, kernelIndCPTV, d_rXZ);
			KARG(kArgs, kernelProxTGVq, kernelIndCPTV, d_rYZ);
		}
		KARG(kArgs, kernelProxTGVq, kernelIndCPTV, alpha);
		// Compute the kernel
#if defined(CUDA) || defined(HIP)
		status = cuLaunchKernel(kernelProxTGVq, globalQ / 64ULL, 1, 1, 64, 1, 1, 0, CLCommandQueue[0], kArgs.data(), NULL);
		CUDA_CHECK(status, "Failed to launch the Proximal TGV kernel\n", -1);
#elif defined(OPENCL)
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelProxTGVq, cl::NullRange, globalQ, cl::NullRange);
		OCL_CHECK(status, "Failed to launch the Proximal TGV kernel\n", -1);
#endif // END CUDA
		//FINISH_QUEUE(status, "Queue finish failed after TGV kernel\n", -1);
		return status;
	}

	/// <summary>
	/// Divergence of the TV prior
	/// </summary>
	/// <param name="im the input array from where the divergence is computed"></param>
	/// <param name="input the backprojection, to which the divergence is added"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="alpha the regularization parameter"></param>
	/// <returns></returns>
#if defined(CUDA) || defined(HIP)
	inline int ProxTVDiv(const scalarStruct & inputScalars) {
#elif defined(OPENCL)
	inline int ProxTVDiv(const scalarStruct & inputScalars, uint32_t timestep = 0) { // TODO: remove default argument
#endif // END CUDA
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrint("Starting Proximal TV divergence");
		STATUS_t status = SUCCESS_VALUE;
#if defined(OPENCL)
		UINT32_t kernelIndCPTV = 0U;
#endif // END CUDA
#if defined(CUDA) || defined(HIP)
		std::vector<void*> kArgs;
#endif // END CUDA
		if (inputScalars.largeDim)
			SET_RANGE_Z(globalPriorEFOV, inputScalars.Nz[0]);
		if (DEBUG) {
			mexPrintBase("erotusPrior[0] = %u\n", erotusPrior[0]);
			mexPrintBase("erotusPrior[1] = %u\n", erotusPrior[1]);
			mexPrintBase("erotusPrior[2] = %u\n", erotusPrior[2]);
			mexPrintBase("globalPriorEFOV[0] = %u\n", globalPriorEFOV[0]);
			mexPrintBase("globalPriorEFOV[1] = %u\n", globalPriorEFOV[1]);
			mexPrintBase("globalPriorEFOV[2] = %u\n", globalPriorEFOV[2]);
			mexPrintBase("d_N.s[0] = %u\n", VEC_X(d_N[0]));
			mexPrintBase("d_N.s[1] = %u\n", VEC_Y(d_N[0]));
			mexPrintBase("d_N.s[2] = %u\n", VEC_Z(d_N[0]));
			mexEval();
		}
		//FINISH_QUEUE(status, "Queue finish failed before divergence kernel\n", -1);
		KARG(kArgs, kernelProxTVDiv, kernelIndCPTV, d_N[0]);
		KARG(kArgs, kernelProxTVDiv, kernelIndCPTV, d_NPrior);
		KARG(kArgs, kernelProxTVDiv, kernelIndCPTV, d_qX);
		KARG(kArgs, kernelProxTVDiv, kernelIndCPTV, d_qY);
		KARG(kArgs, kernelProxTVDiv, kernelIndCPTV, d_qZ);
		KARG(kArgs, kernelProxTVDiv, kernelIndCPTV, vec_opencl.d_rhs_os[0]);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.useBuffers) {
				KARG(kArgs, kernelProxTVDiv, kernelIndCPTV, d_maskPriorB);
			}
			else
#if defined(OPENCL)
				if (inputScalars.maskBPZ > 1) {
					KARG(kArgs, kernelProxTVDiv, kernelIndCPTV, d_maskPrior3);
				}
				else
#endif // END CUDA
					KARG(kArgs, kernelProxTVDiv, kernelIndCPTV, d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			KARG(kArgs, kernelProxTVDiv, kernelIndCPTV, d_eFOVIndices);
		// Compute the kernel
#if defined(CUDA) || defined(HIP)
		status = cuLaunchKernel(kernelProxTVDiv, globalPriorEFOV[0], globalPriorEFOV[1], globalPriorEFOV[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		CUDA_CHECK(status, "Failed to launch the Proximal TV divergence kernel\n", -1);
#elif defined(OPENCL)
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelProxTVDiv, cl::NullRange, globalPriorEFOV, localPrior);
		OCL_CHECK(status, "Failed to launch the Proximal TV divergence kernel\n", -1);
#endif // END CUDA
		//FINISH_QUEUE(status, "Queue finish failed after divergence kernel\n", -1);
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrint("Proximal TV divergence computed");
		return 0;
	}

	/// <summary>
	/// TV prior (gradient)
	/// </summary>
	/// <param name="im the input image (from which the gradient/TV is computed)"></param>
	/// <param name="input the output TV"></param>
	/// <param name="L2Ball if true, computes the projection from an L2 ball, otherwise from the L1 ball"></param>
	/// <param name="sigma adjustable constant for some of the priors"></param>
	/// <param name="v divergence of the symmetric derivative for TGV"></param>
	/// <returns></returns>
#if defined(CUDA) || defined(HIP)
	inline int ProxTVGrad(const scalarStruct & inputScalars, float sigma2, const size_t vSize) {
#elif defined(OPENCL)
	inline int ProxTVGrad(const scalarStruct & inputScalars, const float sigma2, const size_t vSize) {
#endif // END CUDA
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrint("Starting Proximal TV gradient");
		STATUS_t status = SUCCESS_VALUE;
#if defined(OPENCL)
		UINT32_t kernelIndCPTV = 0U;
#endif // END CUDA
#if defined(CUDA) || defined(HIP)
		std::vector<void*> kArgs;
#endif // END CUDA
		if (inputScalars.largeDim)
			SET_RANGE_Z(globalPriorEFOV, inputScalars.Nz[0]);
		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", globalPrior[0]);
			mexPrintBase("global[1] = %u\n", globalPrior[1]);
			mexPrintBase("global[2] = %u\n", globalPrior[2]);
			mexPrintBase("erotusPrior[0] = %u\n", erotusPrior[0]);
			mexPrintBase("erotusPrior[1] = %u\n", erotusPrior[1]);
			mexPrintBase("erotusPrior[2] = %u\n", erotusPrior[2]);
			mexPrintBase("globalPriorEFOV[0] = %u\n", globalPriorEFOV[0]);
			mexPrintBase("globalPriorEFOV[1] = %u\n", globalPriorEFOV[1]);
			mexPrintBase("globalPriorEFOV[2] = %u\n", globalPriorEFOV[2]);
			mexPrintBase("d_N.s[0] = %u\n", VEC_X(d_N[0]));
			mexPrintBase("d_N.s[1] = %u\n", VEC_Y(d_N[0]));
			mexPrintBase("d_N.s[2] = %u\n", VEC_Z(d_N[0]));
			mexPrintBase("vSize = %u\n", vSize);
			mexEval();
		}
		KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, d_N[0]);
		KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, d_NPrior);
		KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, d_inputB);
		KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, d_qX);
		KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, d_qY);
		KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, d_qZ);
		KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, sigma2);
		if (vSize > 0) {
			KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, d_vX);
			KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, d_vY);
			if (!inputScalars.TGV2D)
				KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, d_vZ);
		}
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.useBuffers) {
				KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, d_maskPriorB);
			}
			else
#if defined(OPENCL)
				if (inputScalars.maskBPZ > 1) {
					KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, d_maskPrior3);
				}
				else
#endif // END CUDA
					KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			KARG(kArgs, kernelProxTVGrad, kernelIndCPTV, d_eFOVIndices);
		// Compute the kernel
#if defined(CUDA) || defined(HIP)
		status = cuLaunchKernel(kernelProxTVGrad, globalPriorEFOV[0], globalPriorEFOV[1], globalPriorEFOV[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
#elif defined(OPENCL)
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelProxTVGrad, cl::NullRange, globalPriorEFOV, localPrior);
#endif // END CUDA
		CHECK(status, "Failed to launch the Proximal TV gradient kernel\n", -1);
		if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Proximal TV gradient kernel launched successfully\n");
		}

		//FINISH_QUEUE(status, "Queue finish failed after gradient kernel\n", -1);
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrint("Proximal TV gradient computed");
		return 0;
	}

	/// <summary>
	/// Symmetric derivative for TGV
	/// </summary>
	/// <param name="v input array"></param>
	/// <param name="q the output symmetric derivative array"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="sigma2 the sigma value of CP/PDHG (1 for PKMA)"></param>
	/// <returns></returns>
#if defined(CUDA) || defined(HIP)
	inline int ProxTGVSymmDeriv(const scalarStruct & inputScalars, float sigma2) {
#elif defined(OPENCL)
	inline int ProxTGVSymmDeriv(const scalarStruct & inputScalars, const float sigma2) {
#endif // END CUDA
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrint("Starting Proximal TGV symmetric derivative");
		STATUS_t status = SUCCESS_VALUE;
#if defined(OPENCL)
		UINT32_t kernelIndCPTGV = 0U;
#endif // END CUDA
#if defined(CUDA) || defined(HIP)
		std::vector<void*> kArgs;
		UINT32_t kernelIndCPTGV = 0U;
#endif // END CUDA
		if (inputScalars.largeDim)
			SET_RANGE_Z(globalPriorEFOV, inputScalars.Nz[0]);
		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", globalPriorEFOV[0]);
			mexPrintBase("global[1] = %u\n", globalPriorEFOV[1]);
			mexPrintBase("global[2] = %u\n", globalPriorEFOV[2]);
			mexPrintBase("d_N.s[0] = %u\n", VEC_X(d_N[0]));
			mexPrintBase("d_N.s[1] = %u\n", VEC_Y(d_N[0]));
			mexPrintBase("d_N.s[2] = %u\n", VEC_Z(d_N[0]));
			mexEval();
		}
		KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_N[0]);
		KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_NPrior);
		KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_vX);
		KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_vY);
		if (!inputScalars.TGV2D)
			KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_vZ);
		KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_rX);
		KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_rY);
		if (!inputScalars.TGV2D) {
			KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_rZ);
			KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_rXY);
			KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_rXZ);
			KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_rYZ);
		}
		else
			KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_rXY);
		KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, sigma2);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.useBuffers) {
				KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_maskPriorB);
			}
			else
#if defined(OPENCL)
				if (inputScalars.maskBPZ > 1) {
					KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_maskPrior3);
				}
				else
#endif // END CUDA
					KARG(kArgs, kernelProxTGVSymmDeriv, kernelIndCPTGV, d_maskPrior);
		// Compute the kernel
#if defined(CUDA) || defined(HIP)
		status = cuLaunchKernel(kernelProxTGVSymmDeriv, globalPriorEFOV[0], globalPriorEFOV[1], globalPriorEFOV[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
#elif defined(OPENCL)
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelProxTGVSymmDeriv, cl::NullRange, globalPriorEFOV, localPrior);
#endif // END CUDA
		CHECK(status, "Failed to launch the Proximal TGV symmetric derivative kernel\n", -1);
		if (DEBUG || inputScalars.verbose >= 3) {
			mexPrint("Proximal TGV symmetric derivative kernel launched successfully\n");
		}
		//FINISH_QUEUE(status, "Queue finish failed after symmetric derivative kernel\n", -1);
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrint("Proximal TGV symmetric derivative computed");
		return 0;
	}

	/// <summary>
	/// Divergence for TGV
	/// </summary>
	/// <param name="q first half of the input TGV array"></param>
	/// <param name="q2 second half of the input TGV array"></param>
	/// <param name="v output of the divergence"></param>
	/// <param name="p the TV gradient"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="alpha the regularization parameter"></param>
	/// <param name="theta theta value of CP/PDHG or the momentum parameter for PKMA"></param>
	/// <param name="tau tau value of CP/PDHG (1 for PKMA)"></param>
	/// <returns></returns>
#if defined(CUDA) || defined(HIP)
	inline int ProxTGVDiv(const scalarStruct & inputScalars, float theta, float tau) {
#elif defined(OPENCL)
	inline int ProxTGVDiv(const scalarStruct & inputScalars, const float theta, const float tau) {
#endif // END CUDA
		if (inputScalars.verbose >= 3 || DEBUG) {
			mexPrint("Starting Proximal TGV divergence");
		}
		STATUS_t status = SUCCESS_VALUE;
#if defined(OPENCL)
		UINT32_t kernelIndCPTGV = 0U;
#endif // END CUDA
#if defined(CUDA) || defined(HIP)
		std::vector<void*> kArgs;
#endif // END CUDA
		if (inputScalars.largeDim)
			SET_RANGE_Z(globalPriorEFOV, inputScalars.Nz[0]);
		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", globalPriorEFOV[0]);
			mexPrintBase("global[1] = %u\n", globalPriorEFOV[1]);
			mexPrintBase("global[2] = %u\n", globalPriorEFOV[2]);
			mexPrintBase("d_N.s[0] = %u\n", VEC_X(d_N[0]));
			mexPrintBase("d_N.s[1] = %u\n", VEC_Y(d_N[0]));
			mexPrintBase("d_N.s[2] = %u\n", VEC_Z(d_N[0]));
			mexPrintBase("theta = %f\n", theta);
			mexPrintBase("tau = %f\n", tau);
			mexEval();
		}
		//FINISH_QUEUE(status, "\n", -1);
		KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_N[0]);
		KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_NPrior);
		KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_rX);
		KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_rY);
		if (!inputScalars.TGV2D) {
			KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_rZ);
			KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_rXY);
			KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_rXZ);
			KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_rYZ);
		}
		else
			KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_rXY);
		KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_vX);
		KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_vY);
		if (!inputScalars.TGV2D)
			KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_vZ);
		KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_qX);
		KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_qY);
		KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_qZ);
		KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, theta);
		KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, tau);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.useBuffers) {
				KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_maskPriorB);
			}
			else
#if defined(OPENCL)
				if (inputScalars.maskBPZ > 1) {
					KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_maskPrior3);
				}
				else
#endif // END CUDA
					KARG(kArgs, kernelProxTGVDiv, kernelIndCPTGV, d_maskPrior);
		// Compute the kernel
#if defined(CUDA) || defined(HIP)
		status = cuLaunchKernel(kernelProxTGVDiv, globalPriorEFOV[0], globalPriorEFOV[1], globalPriorEFOV[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		CUDA_CHECK(status, "Failed to launch the Proximal TGV divergence kernel\n", -1);
#elif defined(OPENCL)
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelProxTGVDiv, cl::NullRange, globalPriorEFOV, localPrior);
		OCL_CHECK(status, "Failed to launch the Proximal TGV divergence kernel\n", -1);
#endif // END CUDA
		//FINISH_QUEUE(status, "Queue finish failed after divergence kernel\n", -1);
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrint("Proximal TGV divergence complete");
		return 0;
	}

	/// <summary>
	/// In-place element-wise computations, both multiplication and division supported, for either 1D or 2D arrays
	/// </summary>
	/// <param name="vector input array"></param>
	/// <param name="input input and output array"></param>
	/// <param name="mult if true, performs multiplication, otherwise division"></param>
	/// <param name="D2 if true, assumes 2D case, otherwise 1D"></param>
	/// <returns></returns>
#if defined(CUDA) || defined(HIP)
	inline int elementWiseComp(const bool mult, const uint64_t size[], bool D2 = false) {
		const unsigned int gSize[3] = { static_cast<unsigned int>(size[0]), static_cast<unsigned int>(size[1]), static_cast<unsigned int>(size[2]) };
		std::vector<void*> kArgs;
		if (DEBUG) {
			mexPrintBase("gSize[0] = %u\n", gSize[0]);
			mexPrintBase("gSize[1] = %u\n", gSize[1]);
			mexPrintBase("gSize[2] = %u\n", gSize[2]);
			mexEval();
		}
#elif defined(OPENCL)
	inline int elementWiseComp(const bool mult, const uint64_t size[], const bool D2 = false) {
		cl::NDRange gSize = { static_cast<cl::size_type>(size[0]), static_cast<cl::size_type>(size[1]), static_cast<cl::size_type>(size[2]) };
		UINT32_t kernelIndE = 0U;
#endif // END CUDA
		STATUS_t status = SUCCESS_VALUE;
		UCHAR_t D = static_cast<UCHAR_t>(D2);
		//FINISH_QUEUE(status, "Failed to synchronize before element-wise kernel\n", -1);
		if (mult) {
			KARG(kArgs, kernelElementMultiply, kernelIndE, d_vector);
			KARG(kArgs, kernelElementMultiply, kernelIndE, d_input);
			KARG(kArgs, kernelElementMultiply, kernelIndE, D);
			// Compute the kernel
#if defined(CUDA) || defined(HIP)
			status = cuLaunchKernel(kernelElementMultiply, gSize[0], gSize[1], gSize[2], 1, 1, 1, 0, CLCommandQueue[0], kArgs.data(), NULL);
#elif defined(OPENCL)
			status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelElementMultiply, cl::NullRange, gSize, cl::NullRange);
#endif // END CUDA
		}
		else {
			KARG(kArgs, kernelElementDivision, kernelIndE, d_vector);
			KARG(kArgs, kernelElementDivision, kernelIndE, d_input);
			// Compute the kernel
#if defined(CUDA) || defined(HIP)
			status = cuLaunchKernel(kernelElementDivision, gSize[0], gSize[1], gSize[2], 1, 1, 1, 0, CLCommandQueue[0], kArgs.data(), NULL);
#elif defined(OPENCL)
			status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelElementDivision, cl::NullRange, gSize, cl::NullRange);
#endif // END CUDA
		}
		CHECK(status, "Failed to launch the element-wise kernel\n", -1);
		//FINISH_QUEUE(status, "Queue finish failed after element-wise kernel\n", -1);
		return 0;
	}

	/// <summary>
	/// The gradient of hyperbolic prior
	/// </summary>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="sigma adjustable weighting parameter"></param>
	/// <param name="beta regularization parameter"></param>
	/// <returns></returns>
#if defined(CUDA) || defined(HIP)
	inline int hyperGradient(const scalarStruct & inputScalars, float sigma, const Weighting & w_vec, float beta, const int kk = 0) {
#elif defined(OPENCL)
	inline int hyperGradient(const scalarStruct & inputScalars, const float sigma, const Weighting & w_vec, const float beta, const int kk = 0) {
#endif // END CUDA
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrint("Starting " BACKEND_STR " hyperbolic prior gradient computation");
		STATUS_t status = SUCCESS_VALUE;
#if defined(OPENCL)
		if (inputScalars.largeDim)
			SET_RANGE_Z(globalPrior, inputScalars.Nz[0]);
#endif // END CUDA
		TimerPoint tStart, tEnd;
		if (DEBUG || inputScalars.verbose >= 3) {
			INIT_TIMER(tStart, tEnd);
		}
		if (DEBUG) {
			mexPrintBase("sigma = %f\n", sigma);
			mexPrintBase("beta = %f\n", beta);
			mexEval();
		}
		uint32_t Nz, NzOrig;
		UINT2_t nOffset;
		if (inputScalars.largeDim)
			Nz = inputScalars.lDimStruct.NzPr[kk];
		else
			Nz = inputScalars.Nz[0];
		if (inputScalars.largeDim) {
			SET_RANGE_Z(globalPrior, Nz);
			NzOrig = VEC_Z(d_N[0]);
			VEC_Z(d_N[0]) = Nz;
		}
		if (kk == 0 && inputScalars.largeDim)
			nOffset = { 0, Nz - w_vec.Ndz };
		else if (kk < inputScalars.subsetsUsed - 1 && kk > 0 && inputScalars.largeDim)
			nOffset = { w_vec.Ndz, Nz - w_vec.Ndz };
		else if (inputScalars.largeDim)
			nOffset = { w_vec.Ndz, Nz };
#if defined(CUDA) || defined(HIP)
		std::vector<void*> kArgs;
#elif defined(OPENCL)
		UINT32_t kernelIndHyper = 0U;
#endif // END CUDA
		//FINISH_QUEUE(status, "\n", -1);
		KARG(kArgs, kernelHyper, kernelIndHyper, d_W);
		if (inputScalars.useImages) {
			// FPType == 5 caches integral images rather than the estimate, so the prior uses its own copy
			if (inputScalars.FPType == 5) {
				KARG(kArgs, kernelHyper, kernelIndHyper, d_inputI);
			}
			else {
#if defined(CUDA) || defined(HIP)
				KARG(kArgs, kernelHyper, kernelIndHyper, FPTexCachePrior[kk]);
#else
				KARG(kArgs, kernelHyper, kernelIndHyper, d_imageCache[kk]);
#endif
			}
		}
		else {
			KARG(kArgs, kernelHyper, kernelIndHyper, d_inputB);
		}
#if defined(CUDA) || defined(HIP)
		float smooth = inputScalars.epps;
#endif // END CUDA
		KARG(kArgs, kernelHyper, kernelIndHyper, d_N[0]);
		KARG(kArgs, kernelHyper, kernelIndHyper, d_NOrig);
		KARG(kArgs, kernelHyper, kernelIndHyper, sigma);
#if defined(CUDA) || defined(HIP)
		KARG(kArgs, kernelHyper, kernelIndHyper, smooth);
#elif defined(OPENCL)
		KARG(kArgs, kernelHyper, kernelIndHyper, inputScalars.epps);
#endif // END CUDA
		KARG(kArgs, kernelHyper, kernelIndHyper, beta);
		KARG(kArgs, kernelHyper, kernelIndHyper, d_weights);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.useBuffers) {
				KARG(kArgs, kernelHyper, kernelIndHyper, d_maskPriorB);
			}
			else
#if defined(OPENCL)
				if (inputScalars.maskBPZ > 1) {
					KARG(kArgs, kernelHyper, kernelIndHyper, d_maskPrior3);
				}
				else
#endif // END CUDA
					KARG(kArgs, kernelHyper, kernelIndHyper, d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			KARG(kArgs, kernelHyper, kernelIndHyper, d_eFOVIndices);
		if (inputScalars.largeDim)
			KARG(kArgs, kernelHyper, kernelIndHyper, nOffset);
		// Compute the kernel
		if (DEBUG || inputScalars.verbose >= 3)
			START_TIMER(tStart);
#if defined(CUDA) || defined(HIP)
		status = cuLaunchKernel(kernelHyper, globalPrior[0], globalPrior[1], globalPrior[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		CUDA_CHECK(status, "Failed to launch the hyperbolic prior gradient kernel\n", -1);

		//FINISH_QUEUE(status, "Queue finish failed after hyperbolic prior gradient kernel\n", -1);
		if (inputScalars.useImages && inputScalars.FPType == 5) {
			status = cuTexObjectDestroy(d_inputI);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
			}
			status = cuArrayDestroy(imArray);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
			}
		}
#elif defined(OPENCL)
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelHyper, cl::NullRange, globalPrior, localPrior);
		OCL_CHECK(status, "Failed to launch the hyperbolic prior gradient kernel\n", -1);
		//FINISH_QUEUE(status, "Queue finish failed after hyperbolic prior gradient kernel\n", -1);
#endif // END CUDA
		if (DEBUG || inputScalars.verbose >= 3) {
			STOP_TIMER(tEnd);
			PRINT_TIMER(tStart, tEnd, BACKEND_STR " hyperbolic prior gradient completed in %f seconds\n");
		}
		if (inputScalars.largeDim)
			VEC_Z(d_N[0]) = NzOrig;
		return 0;
	}

	/// <summary>
	/// The gradient of TV prior
	/// </summary>
	/// <param name="grad output gradient array"></param>
	/// <param name="im input image (from which the gradient is computed)"></param>
	/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
	/// <param name="sigma various adjustable parameters for some of the priors"></param>
	/// <param name="smooth smoothing value that allows differentiation"></param>
	/// <returns></returns>
#if defined(CUDA) || defined(HIP)
	inline int TVGradient(const scalarStruct & inputScalars, float sigma, float smooth, const Weighting & w_vec, float beta, const int kk = 0, float C = 0.f, const int type = 0) {
#elif defined(OPENCL)
	inline int TVGradient(const scalarStruct & inputScalars, const float sigma, const float smooth, const Weighting & w_vec, const float beta, const int kk = 0, const float C = 0.f, const int type = 0) {
#endif // END CUDA
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrint("Starting " BACKEND_STR " TV gradient computation");
		STATUS_t status = SUCCESS_VALUE;
#if defined(OPENCL)
		if (inputScalars.largeDim)
			SET_RANGE_Z(globalPrior, inputScalars.Nz[0]);
#endif // END CUDA
		TimerPoint tStart, tEnd;
		if (DEBUG || inputScalars.verbose >= 3) {
			INIT_TIMER(tStart, tEnd);
		}
		//FINISH_QUEUE(status, "\n", -1);
		if (DEBUG) {
			mexPrintBase("sigma = %f\n", sigma);
			mexPrintBase("smooth = %f\n", smooth);
			mexPrintBase("beta = %f\n", beta);
			if (type == 2 || type == 3)
				mexPrintBase("C = %f\n", C);
			mexEval();
		}
		uint32_t Nz, NzOrig;
		UINT2_t nOffset;
		if (inputScalars.largeDim)
			Nz = inputScalars.lDimStruct.NzPr[kk];
		else
			Nz = inputScalars.Nz[0];
		if (inputScalars.largeDim) {
			SET_RANGE_Z(globalPrior, Nz);
			NzOrig = VEC_Z(d_N[0]);
			VEC_Z(d_N[0]) = Nz;
		}
		if (kk == 0 && inputScalars.largeDim)
			nOffset = { 0, Nz - 1 };
		else if (kk < inputScalars.subsetsUsed - 1 && kk > 0 && inputScalars.largeDim)
			nOffset = { 1, Nz - 1 };
		else if (inputScalars.largeDim)
			nOffset = { 1, Nz };
#if defined(CUDA) || defined(HIP)
		std::vector<void*> kArgs;
#elif defined(OPENCL)
		UINT32_t kernelIndTV = 0U;
#endif // END CUDA
		KARG(kArgs, kernelTV, kernelIndTV, d_W);
		if (inputScalars.useImages) {
			// FPType == 5 caches integral images rather than the estimate, so the prior uses its own copy
			if (inputScalars.FPType == 5) {
				KARG(kArgs, kernelTV, kernelIndTV, d_inputI);
			}
			else {
#if defined(CUDA) || defined(HIP)
				KARG(kArgs, kernelTV, kernelIndTV, FPTexCachePrior[kk]);
#else
				KARG(kArgs, kernelTV, kernelIndTV, d_imageCache[kk]);
#endif
			}
		}
		else {
			KARG(kArgs, kernelTV, kernelIndTV, d_inputB);
		}
		KARG(kArgs, kernelTV, kernelIndTV, d_N[0]);
		KARG(kArgs, kernelTV, kernelIndTV, d_NOrig);
		KARG(kArgs, kernelTV, kernelIndTV, sigma);
		KARG(kArgs, kernelTV, kernelIndTV, smooth);
		KARG(kArgs, kernelTV, kernelIndTV, beta);
		if (inputScalars.maskBP || (inputScalars.useExtendedFOV && !inputScalars.multiResolution))
			if (inputScalars.useBuffers) {
				KARG(kArgs, kernelTV, kernelIndTV, d_maskPriorB);
			}
			else
#if defined(OPENCL)
				if (inputScalars.maskBPZ > 1) {
					KARG(kArgs, kernelTV, kernelIndTV, d_maskPrior3);
				}
				else
#endif // END CUDA
					KARG(kArgs, kernelTV, kernelIndTV, d_maskPrior);
		if (inputScalars.eFOV && !inputScalars.multiResolution)
			KARG(kArgs, kernelTV, kernelIndTV, d_eFOVIndices);
		if (type == 2 || type == 3)
			KARG(kArgs, kernelTV, kernelIndTV, C);
		if (type > 0)
			KARG(kArgs, kernelTV, kernelIndTV, d_refIm);
		if (inputScalars.largeDim)
			KARG(kArgs, kernelTV, kernelIndTV, nOffset);
		// Compute the kernel
		if (DEBUG || inputScalars.verbose >= 3)
			START_TIMER(tStart);
#if defined(CUDA) || defined(HIP)
		status = cuLaunchKernel(kernelTV, globalPrior[0], globalPrior[1], globalPrior[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		CUDA_CHECK(status, "Failed to launch the TV gradient kernel\n", -1);
		//FINISH_QUEUE(status, "Queue finish failed after TV gradient kernel\n", -1);
		if (inputScalars.useImages && inputScalars.FPType == 5) {
			status = cuTexObjectDestroy(d_inputI);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
			}
			status = cuArrayDestroy(imArray);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
			}
		}
#elif defined(OPENCL)
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelTV, cl::NullRange, globalPrior, localPrior);
		OCL_CHECK(status, "Failed to launch the TV gradient kernel\n", -1);
		//FINISH_QUEUE(status, "Queue finish failed after TV gradient kernel\n", -1);
#endif // END CUDA
		if (DEBUG || inputScalars.verbose >= 3) {
			STOP_TIMER(tEnd);
			PRINT_TIMER(tStart, tEnd, BACKEND_STR " TV gradient completed in %f seconds\n");
		}
		if (inputScalars.largeDim)
			VEC_Z(d_N[0]) = NzOrig;
		return 0;
	}


#if defined(CUDA) || defined(HIP)
	inline int PoissonUpdate(const scalarStruct & inputScalars, float lambda, float epps, float alpha, const int ii = 0) {
#elif defined(OPENCL)
	inline int PoissonUpdate(const scalarStruct & inputScalars, const float lambda, const float epps, const float alpha, const int ii = 0) {
#endif // END CUDA
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrint("Starting " BACKEND_STR " Poisson update (PKMA/MBSREM/BSREM) computation");
		STATUS_t status = SUCCESS_VALUE;
		TimerPoint tStart, tEnd;
		if (DEBUG || inputScalars.verbose >= 3) {
			INIT_TIMER(tStart, tEnd);
		}
#if defined(CUDA) || defined(HIP)
		std::vector<void*> kArgs;
#elif defined(OPENCL)
		UINT32_t kernelIndPoisson = 0U;
#endif // END CUDA
		//FINISH_QUEUE(status, "\n", -1);
		SET_LAUNCH_RANGE3(global,
			inputScalars.Nx[ii] + erotusPDHG[0][ii],
			inputScalars.Ny[ii] + erotusPDHG[1][ii],
			inputScalars.Nz[ii],
			localPrior);
		UCHAR_t enforcePositivity = static_cast<UCHAR_t>(inputScalars.enforcePositivity);
		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", global[0]);
			mexPrintBase("global[1] = %u\n", global[1]);
			mexPrintBase("global[2] = %u\n", global[2]);
			mexPrintBase("erotusBP[0] = %u\n", erotusBP[0][ii]);
			mexPrintBase("erotusBP[1] = %u\n", erotusBP[1][ii]);
			mexPrintBase("erotusPDHG[0] = %u\n", erotusPDHG[0][ii]);
			mexPrintBase("erotusPDHG[1] = %u\n", erotusPDHG[1][ii]);
			mexPrintBase("localPrior[0] = %u\n", localPrior[0]);
			mexPrintBase("localPrior[1] = %u\n", localPrior[1]);
			mexPrintBase("d_N.s[0] = %u\n", VEC_X(d_N[ii]));
			mexPrintBase("d_N.s[1] = %u\n", VEC_Y(d_N[ii]));
			mexPrintBase("d_N.s[2] = %u\n", VEC_Z(d_N[ii]));
			mexPrintBase("lambda = %.8f\n", lambda);
			mexPrintBase("alpha = %f\n", alpha);
			mexEval();
		}
		KARG(kArgs, kernelPoisson, kernelIndPoisson, d_im);
		KARG(kArgs, kernelPoisson, kernelIndPoisson, d_rhs);
		KARG(kArgs, kernelPoisson, kernelIndPoisson, d_N[ii]);
		KARG(kArgs, kernelPoisson, kernelIndPoisson, lambda);
		KARG(kArgs, kernelPoisson, kernelIndPoisson, epps);
		KARG(kArgs, kernelPoisson, kernelIndPoisson, alpha);
		KARG(kArgs, kernelPoisson, kernelIndPoisson, enforcePositivity);
		// Compute the kernel
		if (DEBUG || inputScalars.verbose >= 3)
			START_TIMER(tStart);
#if defined(CUDA) || defined(HIP)
		status = cuLaunchKernel(kernelPoisson, global[0], global[1], global[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		CUDA_CHECK(status, "Failed to launch the Poisson update kernel\n", -1);
#elif defined(OPENCL)
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelPoisson, cl::NullRange, global, localPrior);
		OCL_CHECK(status, "Failed to launch the Poisson update kernel\n", -1);
#endif // END CUDA
		//FINISH_QUEUE(status, "Queue finish failed after Poisson update kernel\n", -1);
		if (DEBUG || inputScalars.verbose >= 3) {
			STOP_TIMER(tEnd);
			PRINT_TIMER(tStart, tEnd, BACKEND_STR " Poisson update completed in %f seconds\n");
		}
		return 0;
	}

#endif // END non-Metal auxiliary kernels
#if defined(CUDA) || defined(HIP)
	inline int PDHGUpdate(const scalarStruct & inputScalars, float epps, float theta, float tau, const int ii = 0) {
#elif defined(OPENCL) || defined(METAL)
	inline int PDHGUpdate(const scalarStruct & inputScalars, const float epps, const float theta, const float tau, const int ii = 0) {
#endif // END CUDA
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrint("Starting " BACKEND_STR " PDHG update computation");
		STATUS_t status = SUCCESS_VALUE;
		TimerPoint tStart, tEnd;
		if (DEBUG || inputScalars.verbose >= 3) {
			INIT_TIMER(tStart, tEnd);
		}
#if defined(CUDA) || defined(HIP)
		std::vector<void*> kArgs;
#elif defined(OPENCL) || defined(METAL)
		UINT32_t kernelIndPDHG = 0U;
#endif // END CUDA
#if defined(METAL)
		if (!queueBP || !kernelPDHG) {
			mexPrint("Unable to create Metal PDHG update encoder");
			return -1;
		}
		NS::SharedPtr<MTL::CommandBuffer> commandBuffer = NS::RetainPtr(queueBP->commandBuffer());
		NS::SharedPtr<MTL::ComputeCommandEncoder> encoder = NS::RetainPtr(commandBuffer->computeCommandEncoder());
		if (!commandBuffer || !encoder) {
			mexPrint("Unable to create Metal PDHG update encoder");
			return -1;
		}
		encoder->setComputePipelineState(kernelPDHG.get());
#endif
		FINISH_QUEUE(status, "\n", -1);
		SET_LAUNCH_RANGE3(global,
			inputScalars.Nx[ii] + erotusPDHG[0][ii],
			inputScalars.Ny[ii] + erotusPDHG[1][ii],
			inputScalars.Nz[ii],
			localPrior);
		UCHAR_t enforcePositivity = static_cast<UCHAR_t>(inputScalars.enforcePositivity);
		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", global[0]);
			mexPrintBase("global[1] = %u\n", global[1]);
			mexPrintBase("global[2] = %u\n", global[2]);
			mexPrintBase("d_N.s[0] = %u\n", VEC_X(d_N[ii]));
			mexPrintBase("d_N.s[1] = %u\n", VEC_Y(d_N[ii]));
			mexPrintBase("d_N.s[2] = %u\n", VEC_Z(d_N[ii]));
			mexPrintBase("theta = %f\n", theta);
			mexPrintBase("tau = %f\n", tau);
			mexEval();
		}
		KARG(kArgs, kernelPDHG, kernelIndPDHG, d_im);
		KARG(kArgs, kernelPDHG, kernelIndPDHG, d_rhs);
		KARG(kArgs, kernelPDHG, kernelIndPDHG, d_U);
#if defined(METAL)
		kParams.N_PDHG = { static_cast<int>(VEC_X(d_N[ii])), static_cast<int>(VEC_Y(d_N[ii])), static_cast<int>(VEC_Z(d_N[ii])) };
		kParams.epps_PDHG = epps;
		kParams.theta_PDHG = theta;
		kParams.tau_PDHG = tau;
		kParams.enforcePositivity_PDHG = enforcePositivity;
		KARG(kArgs, kernelPDHG, kernelIndPDHG, kParams);
#else
		KARG(kArgs, kernelPDHG, kernelIndPDHG, d_N[ii]);
		KARG(kArgs, kernelPDHG, kernelIndPDHG, epps);
		KARG(kArgs, kernelPDHG, kernelIndPDHG, theta);
		KARG(kArgs, kernelPDHG, kernelIndPDHG, tau);
		KARG(kArgs, kernelPDHG, kernelIndPDHG, enforcePositivity);
#endif
		// Compute the kernel
		if (DEBUG || inputScalars.verbose >= 3)
			START_TIMER(tStart);
#if defined(CUDA) || defined(HIP)
		status = cuLaunchKernel(kernelPDHG, global[0], global[1], global[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		CUDA_CHECK(status, "Failed to launch the PDHG update kernel\n", -1);
#elif defined(METAL)
		{
			const MTL::Size threadsPerThreadgroup = MTL::Size::Make(localPrior[0], localPrior[1], localPrior[2]);
			const MTL::Size threadgroupsPerGrid = MTL::Size::Make(
				global[0] / localPrior[0], global[1] / localPrior[1], global[2] / localPrior[2]);
			encoder->dispatchThreadgroups(threadgroupsPerGrid, threadsPerThreadgroup);
			encoder->endEncoding();
			commandBuffer->commit();
		}
#elif defined(OPENCL)
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelPDHG, cl::NullRange, global, localPrior);
		OCL_CHECK(status, "Failed to launch the PDHG update kernel\n", -1);
#endif // END CUDA
		//FINISH_QUEUE(status, "Queue finish failed after PDHG update kernel\n", -1);
		if (DEBUG || inputScalars.verbose >= 3) {
			STOP_TIMER(tEnd);
			PRINT_TIMER(tStart, tEnd, BACKEND_STR " PDHG update completed in %f seconds\n");
		}
		return 0;
	}

#if defined(CUDA) || defined(HIP)
	inline int rotateCustom(const scalarStruct & inputScalars, float cosa, float sina, const int ii = 0) {
#elif defined(OPENCL) || defined(METAL)
	inline int rotateCustom(const scalarStruct & inputScalars, const float cosa, const float sina, const int ii = 0) {
#endif // END CUDA
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrint("Starting " BACKEND_STR " bilinear image rotation computation");
		STATUS_t status = SUCCESS_VALUE;
#if defined(CUDA) || defined(HIP)
		std::vector<void*> kArgs;
#elif defined(OPENCL) || defined(METAL)
		UINT32_t kernelIndRot = 0U;
#endif // END CUDA
#if defined(METAL)
		if (!queueBP || !kernelRotate) {
			mexPrint("Unable to create Metal bilinear rotation encoder");
			return -1;
		}
		NS::SharedPtr<MTL::CommandBuffer> commandBuffer = NS::RetainPtr(queueBP->commandBuffer());
		if (!commandBuffer) {
			mexPrint("Unable to create Metal bilinear rotation encoder");
			return -1;
		}
		NS::SharedPtr<MTL::ComputeCommandEncoder> encoder = NS::RetainPtr(commandBuffer->computeCommandEncoder());
		if (!encoder) {
			mexPrint("Unable to create Metal bilinear rotation encoder");
			return -1;
		}
		encoder->setComputePipelineState(kernelRotate.get());
#endif
		SET_LAUNCH_RANGE3(global,
			inputScalars.Nx[ii] + erotusPrior[0],
			inputScalars.Ny[ii] + erotusPrior[1],
			inputScalars.Nz[ii],
			localPrior);
		if (DEBUG) {
			mexPrintBase("global[0] = %u\n", global[0]);
			mexPrintBase("global[1] = %u\n", global[1]);
			mexPrintBase("global[2] = %u\n", global[2]);
			mexPrintBase("d_N.s[0] = %u\n", VEC_X(d_N[ii]));
			mexPrintBase("d_N.s[1] = %u\n", VEC_Y(d_N[ii]));
			mexPrintBase("d_N.s[2] = %u\n", VEC_Z(d_N[ii]));
			mexEval();
		}

		KARG(kArgs, kernelRotate, kernelIndRot, d_rhs);
		if (!inputScalars.useBuffers) {
			KARG(kArgs, kernelRotate, kernelIndRot, d_inputI);
		}
		else {
			KARG(kArgs, kernelRotate, kernelIndRot, d_im);
		}
#if defined(METAL)
		kParams.N_rotate = { static_cast<int>(VEC_X(d_N[ii])), static_cast<int>(VEC_Y(d_N[ii])), static_cast<int>(VEC_Z(d_N[ii])) };
		kParams.cosa_rotate = cosa;
		kParams.sina_rotate = sina;
		KARG(kArgs, kernelRotate, kernelIndRot, kParams);
#else
		KARG(kArgs, kernelRotate, kernelIndRot, VEC_X(d_N[ii]));
		KARG(kArgs, kernelRotate, kernelIndRot, VEC_Y(d_N[ii]));
		KARG(kArgs, kernelRotate, kernelIndRot, VEC_Z(d_N[ii]));
		KARG(kArgs, kernelRotate, kernelIndRot, cosa);
		KARG(kArgs, kernelRotate, kernelIndRot, sina);
#endif
		// Compute the kernel
#if defined(CUDA) || defined(HIP)
		status = cuLaunchKernel(kernelRotate, global[0], global[1], global[2], localPrior[0], localPrior[1], localPrior[2], 0, CLCommandQueue[0], kArgs.data(), NULL);
		CUDA_CHECK(status, "Failed to launch the bilinear image rotation kernel\n", -1);
		//FINISH_QUEUE(status, "Queue finish failed after bilinear image rotation kernel\n", -1);
		if (inputScalars.useImages) {
			status = cuTexObjectDestroy(d_inputI);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
			}
			status = cuArrayDestroy(imArray);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
			}
		}
#elif defined(METAL)
		{
			const MTL::Size threadsPerThreadgroup = MTL::Size::Make(localPrior[0], localPrior[1], localPrior[2]);
			const MTL::Size threadgroupsPerGrid = MTL::Size::Make(
				global[0] / localPrior[0], global[1] / localPrior[1], global[2] / localPrior[2]);
			encoder->dispatchThreadgroups(threadgroupsPerGrid, threadsPerThreadgroup);
			encoder->endEncoding();
			commandBuffer->commit();
		}
#elif defined(OPENCL)
		status = (CLCommandQueue[0]).enqueueNDRangeKernel(kernelRotate, cl::NullRange, globalPrior, localPrior);
		OCL_CHECK(status, "Failed to launch the bilinear image rotation kernel\n", -1);
		//FINISH_QUEUE(status, "Queue finish failed after bilinear image rotation kernel\n", -1);
#endif // END CUDA
		if (inputScalars.verbose >= 3)
			mexPrint(BACKEND_STR " bilinear image rotation computed");
		return 0;
		}

#if defined(CUDA) || defined(HIP) || defined(METAL)
	inline int transferTex(const scalarStruct& inputScalars, AFDEVBUFF_t input, const bool RDP = false, const uint32_t Nz = 1) {
		STATUS_t status = SUCCESS_VALUE;
		if (RDP)
			CREATE_FLOAT_TEXTURE3D_FROM_DEVICE(d_RDPrefI, imArray, input, inputScalars.Nx[0], inputScalars.Ny[0], Nz,
				BACKEND_TEXTURE_POINT, BACKEND_TEXTURE_DEFAULT_FLAGS);
		else
			CREATE_FLOAT_TEXTURE3D_FROM_DEVICE(d_inputI, imArray, input, inputScalars.Nx[0], inputScalars.Ny[0], Nz,
				BACKEND_TEXTURE_POINT, BACKEND_TEXTURE_DEFAULT_FLAGS);
		CHECK(status, "Image copy failed\n", -1);
		FINISH_QUEUE(status, "Synchronization failed\n", -1);
		if (DEBUG)
			mexPrint("Synchronization completed\n");
		return 0;
	}
#endif // END CUDA/METAL
	};
