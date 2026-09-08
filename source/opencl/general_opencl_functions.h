/*******************************************************************************************************************************************
* General functions for all the OpenCL and CUDA kernel files. Contains functions that compute the necessary source and detector coordinates, 
* atomics, forward, backward projections, etc. Special functions are available for different cases such as TOF, listmode data, CT data, 
* etc.
*
* Note that all functions are used for both OpenCL and CUDA. To achieve this, preprocessor definitions are used VERY extensively. This can
* make following the code sometimes difficult. The start of the file contains the preprocessor definitions for OpenCL and then for CUDA.
* Note that these definitions are also used in the projector kernel files and in the "auxliary" kernel file.
*
* USEIMAGES specifies whether OpenCL images or CUDA textures are used. If it is not defined, regular buffers are used. Default is ON.
*
* Copyright (C) 2019-2026 Ville-Veikko Wettenhovi, Niilo Saarlemo
*
* This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************************************************************************/

// hipRTC force-includes hiprtc_runtime.h, which uses FP as a type name; defining FP on the
// command line breaks that header. HIP RTC builds therefore pass -DOMEGA_FP instead, and it
// is mapped back to FP here, after the hipRTC prelude has already been processed.
#if defined(OMEGA_FP) && !defined(FP)
#define FP
#endif
#ifdef ATOMIC
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable
#endif
#ifndef N_REKOS
#define N_REKOS 1
#endif
#define NROLLS (N_REKOS * NBINS)
#ifdef PRECOMPUTE
#define TYPE 1
#else
#define TYPE 0
#endif
#ifdef VOL
#define CC 1e3f
#endif
#define TRAPZ_BINS 4.f
#ifdef PITCH
#define NA 6
#else
#define NA 2
#endif
#if defined(PTYPE4) && !defined(BP4)
#define typeT float3
#define T4 float4
#define typeTT float
#else
#ifdef USEIMAGES
#define typeT int3
#else
#if defined(OPENCL) || defined(METAL)
#define typeT long
#else
#define typeT long long
#endif
#endif
#define T4 int4
#define typeTT int
#endif
// The backprojection mask is read with normalized float coordinates only in the ray-based
// (non-CT) projector type 4 backprojection
#if defined(MASKBP) && defined(PTYPE4) && !defined(BP4) && !defined(CT) && (defined(BP) || defined(SENS))
#define MASKBPNORM
#define T3 float3
#define T2 float2
#else
#define T3 int3
#define T2 int2
#endif
// FP mask is always read with integer indices
#ifdef USEIMAGES
#define MASKFPT int3
#else
#if defined(OPENCL) || defined(METAL)
#define MASKFPT long
#else
#define MASKFPT long long
#endif
#endif

#ifdef HALF
#define FLOAT_ZERO 0.h
#define FLOAT_HALF .5h
#define FLOAT_ONE 1.h
#define FLOAT_TWO 2.h
#define FLOAT half
#define FLOAT2 half2
#define FLOAT3 half3
#else
#define FLOAT_ZERO 0.f
#define FLOAT_HALF .5f
#define FLOAT_ONE 1.f
#define FLOAT_TWO 2.f
#define FLOAT float
#define FLOAT2 float2
#define FLOAT3 float3
#endif

// Macro to unpack scalar parameters for projector types 1, 2 and 3 (FP+BP)
#define UNPACK_SCALAR_PARAMS_123(scalarParams) \
    FLOAT global_factor = scalarParams.global_factor; \
    FLOAT d_epps = scalarParams.epps; \
	uint d_size_x = scalarParams.nRowsD; \
	uint d_det_per_ring = scalarParams.det_per_ring; \
	FLOAT sigma_x = scalarParams.sigma_x; \
	FLOAT coneOfResponseStdCoeffA = scalarParams.coneOfResponseStdCoeffA; \
    FLOAT coneOfResponseStdCoeffB = scalarParams.coneOfResponseStdCoeffB; \
    FLOAT coneOfResponseStdCoeffC = scalarParams.coneOfResponseStdCoeffC; \
	FLOAT2 crystalSize = scalarParams.dPitch; \
    FLOAT3 ellipseCenter = scalarParams.ellipseCenter; \
    FLOAT3 ellipseRadii = scalarParams.ellipseRadii; \
    FLOAT ellipsePower = scalarParams.ellipsePower; \
	FLOAT bmin = scalarParams.bmin; \
	FLOAT bmax = scalarParams.bmax; \
	FLOAT Vmax = scalarParams.Vmax; \
	uint d_sizey = scalarParams.nColsD; \
	long d_nProjections = scalarParams.nProjections; \
	uint rings = scalarParams.rings; \
	uint3 d_Nxyz = scalarParams.d_N; \
	FLOAT3 d_d = scalarParams.d; \
	FLOAT3 b = scalarParams.b; \
	FLOAT3 d_bmax = scalarParams.d_bmax; \
	unsigned char no_norm = scalarParams.no_norm; \
	unsigned long m_size = scalarParams.m_size; \
	uint currentSubset = scalarParams.currentSubset; \
	int aa = scalarParams.aa; \
    float orthWidth = scalarParams.orthWidth;

#define UNPACK_SCALAR_PARAMS_4_FP(scalarParams) \
    const uint d_size_x = scalarParams.nRowsD; \
    const uint d_sizey = scalarParams.nColsD; \
    const float2 d_dPitch = scalarParams.dPitch; \
    const float helicalRadius = scalarParams.helicalRadius; \
    const float dL = scalarParams.dL; \
    const float global_factor = scalarParams.global_factor; \
    const float sigma_x = scalarParams.sigma_x; \
    const uint3 d_N = scalarParams.d_N; \
    const float3 b = scalarParams.b; \
    const float3 bmax = scalarParams.d_bmax; \
    const float3 d_scale = scalarParams.d_Scale4; \
    const int rings = scalarParams.rings; \
    const uint d_det_per_ring = scalarParams.det_per_ring; \
    const long d_nProjections = scalarParams.nProjections; \
    const uchar no_norm = scalarParams.no_norm; \
    const unsigned long m_size = scalarParams.m_size; \
    const uint currentSubset = scalarParams.currentSubset; \
	const int aa = scalarParams.aa;

#define UNPACK_SCALAR_PARAMS_4_BP(scalarParams) \
    const uint d_size_x = scalarParams.nRowsD; \
    const uint d_sizey = scalarParams.nColsD; \
    const float2 d_dPitch = scalarParams.dPitch; \
    const float helicalRadius = scalarParams.helicalRadius; \
    const uint3 d_N = scalarParams.d_N; \
    const float3 b = scalarParams.b; \
    const float3 d_d = scalarParams.d; \
    const float kerroin = scalarParams.kerroin4; \
    const float DSC = scalarParams.DSC; \
    const uchar no_norm = scalarParams.no_norm; \
    const long d_nProjections = scalarParams.nProjections; \
    const int ii = scalarParams.aa;

#define UNPACK_SCALAR_PARAMS_5_FP(scalarParams) \
    const uint d_nRows = scalarParams.nRowsD; \
    const uint d_nCols = scalarParams.nColsD; \
    const float2 d_dPitch = scalarParams.dPitch; \
    const uint3 d_N = scalarParams.d_N; \
    const float3 b = scalarParams.b; \
    const float2 d_Size = scalarParams.dSize5; \
    const float3 d_d = scalarParams.d; \
    const float3 d_scale = scalarParams.d_Scale5; \
    const long d_nProjections = scalarParams.nProjections;

#define UNPACK_SCALAR_PARAMS_PDHG(scalarParams) \
    const int3 N = scalarParams.N_PDHG; \
    const float epps = scalarParams.epps_PDHG; \
    const float theta = scalarParams.theta_PDHG; \
    const float tau = scalarParams.tau_PDHG; \
    const uchar enforcePositivity = scalarParams.enforcePositivity_PDHG;

#define UNPACK_SCALAR_PARAMS_ROTATE(scalarParams) \
    const int3 N_rotate = scalarParams.N_rotate; \
    const int Nx = N_rotate.x; \
    const int Ny = N_rotate.y; \
    const int Nz = N_rotate.z; \
    const float cosa = scalarParams.cosa_rotate; \
    const float sina = scalarParams.sina_rotate;


#ifdef METAL

#if defined(ATOMIC32)
#define CAST int
#ifndef TH
#define TH 100000.f
#endif
#else
#define CAST float
#endif

constexpr metal::sampler samplerForw(
    metal::coord::normalized,
    metal::filter::linear,
    metal::address::clamp_to_edge
);

constexpr metal::sampler sampler2(
    metal::coord::normalized,
    metal::filter::linear,
    metal::address::clamp_to_edge
);

constexpr metal::sampler samplerMask(
    metal::coord::normalized,
    metal::filter::nearest,
    metal::address::clamp_to_edge
);

#ifdef HALF // 16-bit floating point
#define CFLOAT(a) static_cast<half>(a)
#define CFLOAT3(a) half3(a)
#define CMFLOAT3(x,y,z) half3((x), (y), (z))

#define make_float2(a,b) half2((a),(b))
#define make_float3(a,b,c) half3((a),(b),(c)) // TODO: replace with half
#define MFLOAT2(a,b) half2((a), (b)) // TODO: replace with half
#else // 32-bit floating point
#define CFLOAT(a) static_cast<float>(a)
#define CFLOAT3(a) float3(a)
#define CMFLOAT3(x,y,z) float3((x), (y), (z))
#define FLOAT float
#define make_float2(a,b) float2((a),(b))
#define make_float3(a,b,c) float3((a),(b),(c))
#define MFLOAT2(a,b) float2((a), (b))
#define MFLOAT3(a,b,c) float3((a), (b), (c))
#endif

#define ACOS metal::acos
#define ALL metal::all
#define ANY metal::any
#define BUF0 [[buffer(0)]]
#define BUF1 [[buffer(1)]]
#define BUF2 [[buffer(2)]]
#define BUF3 [[buffer(3)]]
#define BUF4 [[buffer(4)]]
#define BUF5 [[buffer(5)]]
#define BUF6 [[buffer(6)]]
#define BUF7 [[buffer(7)]]
#define BUF8 [[buffer(8)]]
#define BUF9 [[buffer(9)]]
#define BUF10 [[buffer(10)]]
#define BUF11 [[buffer(11)]]
#define BUF12 [[buffer(12)]]
#define BUF13 [[buffer(13)]]
#define BUF14 [[buffer(14)]]
#define BUF15 [[buffer(15)]]
#define CEIL metal::ceil
#define CINT(a)   static_cast<int>(a)
#define CINT_rtz(a) static_cast<int>(metal::trunc((a)))
#define CINT3_rtz(a) static_cast<int3>(a)
#define CLAMP3(a,b,c) metal::clamp((a),(b),(c))
#define CLGLOBAL device
#define CLONG_rtz(a) static_cast<long>(a)
#define CLRESTRICT
#define CMINT3(a,b,c) int3((a),(b),(c))
#define CONSTANT constant
#define CROSS metal::cross
#define CUINT(a) (uint)(a)
#define CUINT3(a) uint3(a)
#define CUINT_rtp(a) static_cast<uint>(metal::ceil((a)))
#define CUINT_rtz(a) static_cast<uint>(metal::trunc((a)))
#define CUINT_sat_rtz(a) static_cast<uint>(metal::clamp(metal::trunc(((float)a)), 0.0f, 4294967295.0f)) // TODO replace float with FLOAT
#define DEVICE inline
#define DISTANCE metal::distance
#define DIVIDE(a,b) ((a) / (b))
#define DIVIDE3(a,b) ((a) / (b))
#define EXP(a) metal::exp(a)
#define FABS metal::fabs
#define FMAD(a,b,c) metal::fma(a,b,c)
#define FMAD3(a,b,c) metal::fma((a),(b),(c))
#define FMAX metal::fmax
#define FMIN metal::fmin
#define FLOOR metal::floor
#define IMAGE2D metal::texture2d<float, metal::access::sample>
#define IMAGE3D metal::texture3d<float, metal::access::sample>
#define ISINF metal::isinf
#define ISNAN metal::isnan
#define KERNEL kernel
#define KERNEL2 kernel
#define KERNEL3 kernel
#define LENGTH metal::length
#define LOCAL threadgroup
#define LOG(a) metal::log(a)
#define LONG long
#define make_int2(a,b) int2((a),(b))
#define make_int3(a,b,c) int3((a),(b),(c))
#define make_uint2(a,b) uint2((a),(b))
#define make_uint3(a,b,c) uint3((a),(b),(c))
#define MAX metal::max
#define MIN metal::min
#ifndef LTYPE3
#define MINT3(a,b,c) int3((a),(b),(c))
#else
#define MINT3(a,b,c) long3((a),(b),(c))
#endif
#define MUINT3(a, b, c) uint3(a, b, c)
#define NORMALIZE metal::normalize
#define POWR metal::pow
#define SINF(a) metal::sin(a)
#define COSF(a) metal::cos(a)
#define PTR_DEV device
#define PTR_THR thread
#define PTR_CONST constant
#define PTR_TG threadgroup
#define RCP(x) (1.f / x)
#define SQRT metal::sqrt
#define RSQRT(x) metal::rsqrt(x)
#define SCALAR_PARAMS(name) constant ScalarKernelParams &name
#ifdef USEIMAGES
#define TEX1 [[texture(1)]]
#define TEX2 [[texture(2)]]
#define TEX3 [[texture(3)]]
#define TEX4 [[texture(4)]]
#define TEX7 [[texture(7)]]
#define TEX8 [[texture(8)]]
#define TEX9 [[texture(9)]]
#define TEX19 [[texture(19)]]
#else
#define TEX1 [[buffer(1)]]
#define TEX2 [[buffer(2)]]
#define TEX3 [[buffer(3)]]
#define TEX4 [[buffer(4)]]
#define TEX7 [[buffer(7)]]
#define TEX8 [[buffer(8)]]
#define TEX9 [[buffer(9)]]
#define TEX19 [[buffer(19)]]
#endif
// Metal function definitions
using metal::dot;

#if defined(ATOMIC32)
inline void atomicAdd(volatile device metal::atomic_int* addr, int val)
{
    atomic_fetch_add_explicit(addr, val, metal::memory_order_relaxed);
}
#else
inline void atomicAdd(volatile device metal::atomic_float* addr, float val)
{
    atomic_fetch_add_explicit(addr, val, metal::memory_order_relaxed);
}
#endif

#endif
#ifdef OPENCL
#define BUF0
#define BUF1
#define BUF2
#define BUF3
#define BUF4
#define BUF5
#define BUF6
#define BUF7
#define BUF8
#define BUF9
#define BUF10
#define BUF11
#define BUF12
#define BUF13
#define BUF14
#define BUF15
#define BUF16
#define BUF17
#define BUF18
#define BUF19
#define BUF20
#define TEX0
#define TEX1
#define TEX2
#define TEX3
#define TEX4
#define TEX5
#define TEX6
#define TEX7
#define TEX8
#define TEX9
#define TEX10
#define TEX11
#define TEX12
#define TEX13
#define TEX14
#define TEX15
#define TEX16
#define TEX17
#define TEX18
#define TEX19
#define TEX20
#define PTR_DEV // Metal requires address space qualifier for pointers
#define PTR_THR 
#define PTR_CONST
#define PTR_TG
#define MIN min
#define FABS fabs
#define FMIN fmin
#define FMAX fmax
#define MIN min
#define MAX max
#define ALL all
#define ANY any
#define ISNAN isnan
#define ISINF isinf
#define LONG long
#define ULONG ulong
#define CLGLOBAL __global
#define CLRESTRICT restrict
#define CONSTANT __constant
#define LENGTH length
#define LOCAL __local
#define IMAGE3D __read_only image3d_t
#define IMAGE2D __read_only image2d_t
#define DEVICE
#define UINT_sat(a) convert_uint_sat(a)
#define RIMAGEF(a, b, c) read_imagef(a, b, c).w
#define CINT_rtz(a) convert_int_rtz(a)
#define DIVIDE(a,b) native_divide(a,b)
#define DIVIDE3(a,b) native_divide(a,b)
#define CFLOAT(a) convert_float(a)
#define CFLOAT3(a) convert_float3(a)
#define CUINT(a) convert_uint(a)
#define CUINT3(a) convert_uint3(a)
#define CUINT_rtp(a) convert_uint_rtp(a)
#define CUINT_rtz(a) convert_uint_rtz(a)
#define CUINT_rte(a) convert_uint_rte(a)
#define CUINT_sat_rtz(a) convert_uint_sat_rtz(a)
#define CLONG_rtz(a) convert_long_rtz(a)
#define EXP(a) native_exp(a)
#define EXP3(a) native_exp(a)
#define SINF(a) native_sin(a)
#define COSF(a) native_cos(a)
#define SQRT native_sqrt
#define RSQRT(x) native_rsqrt(x)
#define POWR native_powr
#define POWN pown
#define RCP(x) native_recip(x)
#define FLOOR floor
#define CEIL ceil
#define ATAN2 atan2
#define NORMALIZE normalize
#define CROSS cross
#define DISTANCE distance
#define ACOS acos
#define LOG native_log
#define CLAMP3(a, b, c) clamp(a, b, c)
#define CINT(a) convert_int(a)
#define CINT3(a) convert_int3(a)
#define CINT3_rtz(a) convert_int3_rtz(a)
#define FMAD(a,b,c) mad(a,b,c)
#define FMAD2(a,b,c) mad(a,b,c)
#define FMAD3(a,b,c) mad(a,b,c)
#define GID0 get_global_id(0)
#define GID1 get_global_id(1)
#define GID2 get_global_id(2)
#define GSIZE0 get_global_size(0)
#define GSIZE1 get_global_size(1)
#define GSIZE2 get_global_size(2)
#define GRID0 get_group_id(0)
#define GRID1 get_group_id(1)
#define GRID2 get_group_id(2)
#define LSIZE0 get_local_size(0)
#define LSIZE1 get_local_size(1)
#define LSIZE2 get_local_size(2)
#define LID0 get_local_id(0)
#define LID1 get_local_id(1)
#define LID2 get_local_id(2)
#define CFLOAT4 (float4)
#define CFLOAT2 (float2)
#define MUINT2(a, b) {a, b}
#define MINT3(a, b, c) {a, b, c}
#define MUINT3(a, b, c) {a, b, c}
#define MFLOAT2(a, b) {a, b}
#define MFLOAT3(a, b, c) {a, b, c}
#define CMFLOAT3 (float3)
#define CMINT3 (int3)
#define CMINT4 (int4)
#define BARRIER barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
#define KERNEL __kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, 1)))
#define KERNEL2 __kernel __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, 1)))
#define KERNEL3 __kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, LOCAL_SIZE3)))
#define KERN __kernel
__constant sampler_t samplerIm = CLK_NORMALIZED_COORDS_TRUE | CLK_FILTER_LINEAR | CLK_ADDRESS_CLAMP_TO_EDGE;
#ifdef PTYPE4
__constant sampler_t samplerForw = CLK_NORMALIZED_COORDS_TRUE | CLK_FILTER_LINEAR | CLK_ADDRESS_CLAMP_TO_EDGE;
#endif
__constant sampler_t samplerSiddon = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_CLAMP_TO_EDGE;

#ifdef MASKBPNORM
__constant sampler_t sampler_MASK = CLK_NORMALIZED_COORDS_TRUE | CLK_FILTER_NEAREST | CLK_ADDRESS_CLAMP_TO_EDGE;
#else
__constant sampler_t sampler_MASK = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_CLAMP_TO_EDGE;
#endif
// The forward projection mask is always read with (unnormalized) integer coordinates
__constant sampler_t sampler_MASKFP = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_CLAMP_TO_EDGE;
#endif
// CUDA and HIP share the same device-side language (kernel/qualifier keywords, intrinsics, vector
// types, texture fetch templates, etc.), so the vast majority of these definitions are common. The
// HIP-specific symbols (e.g. the texture object type) are separated where they differ from CUDA.
#if defined(CUDA) || defined(HIP)
#define BUF0
#define BUF1
#define BUF2
#define BUF3
#define BUF4
#define BUF5
#define BUF6
#define BUF7
#define BUF8
#define BUF9
#define BUF10
#define BUF11
#define BUF12
#define BUF13
#define BUF14
#define BUF15
#define BUF16
#define BUF17
#define BUF18
#define BUF19
#define BUF20
#define TEX0
#define TEX1
#define TEX2
#define TEX3
#define TEX4
#define TEX5
#define TEX6
#define TEX7
#define TEX8
#define TEX9
#define TEX10
#define TEX11
#define TEX12
#define TEX13
#define TEX14
#define TEX15
#define TEX16
#define TEX17
#define TEX18
#define TEX19
#define TEX20
#define PTR_DEV // Metal requires address space qualifier for pointers
#define PTR_THR 
#define PTR_CONST
#define PTR_TG
#define MIN min
#define FABS fabs
#define LENGTH length
#define NORMALIZE normalize
#define FMIN fmin
#define FMAX fmax
#define MIN min
#define MAX max
#define ALL all
#define ANY any
#define ISNAN isnan
#define ISINF isinf
#define M_PI_F 3.141593f
#define M_PI_2_F 1.570796f
#define M_SQRT1_2_F 0.7071068f
#define M_1_PI_F 0.3183099f
#define CAST float
#define LONG long long
#define ULONG unsigned long long
#define uint unsigned int
#define ushort unsigned short
#define uchar unsigned char
#define CUINT(a) (unsigned int)(a)
#define UINT_sat(a) (unsigned int)(a)
#define CINT_rtz(a) __float2int_rz(a)
#define DIVIDE(a,b) fdividef(a,b)
#define DIVIDE3(a,b) fdividef3(a,b)
#define CFLOAT(a) (float)(a)
#define CFLOAT3(a) convert_float3(a)
#define CUINT(a) (unsigned int)(a)
#define CUINT3(a) convert_uint3(a)
#define CUINT_rtp(a) __float2uint_ru(a)
#define CUINT_rtz(a) __float2uint_rd(a)
#define CUINT_rte(a) __float2uint_rn(a)
#define CUINT_sat_rtz(a) __float2uint_rd(a)
#define CINT(a) (int)(a)
#define CINT3(a) make_long3_int3(a)
#define CINT3_rtz(a) __float2int_rz3(a)
#define CLONG_rtz(a) __float2ll_rz(a)
#define SINF(a) sinf(a)
#define COSF(a) cosf(a)
#define CLGLOBAL
#define CLRESTRICT
#define CONSTANT const
#define DEVICE inline __device__
#define LOCAL __shared__
#define EXP(a) __expf(a)
#define EXP3(a) expf3(a)
#define FMAD(a, b, c) __fmaf_rn(a, b, c)
#define FMAD2(a, b, c) __fmaf_rn2(a, b, c)
#define FMAD3(a, b, c) __fmaf_rn3(a, b, c)
#define CLAMP3(a, b, c) clamp3(a, b, c)
#define GID0 (threadIdx.x + blockIdx.x * blockDim.x)
#define GID1 (threadIdx.y + blockIdx.y * blockDim.y)
#define GID2 (threadIdx.z + blockIdx.z * blockDim.z)
#define GSIZE0 (blockDim.x * gridDim.x)
#define GSIZE1 (blockDim.y * gridDim.y)
#define GSIZE2 (blockDim.z * gridDim.z)
#define GRID0 blockIdx.x
#define GRID1 blockIdx.y
#define GRID2 blockIdx.z
#define LSIZE0 blockDim.x
#define LSIZE1 blockDim.y
#define LSIZE2 blockDim.z
#define LID0 threadIdx.x
#define LID1 threadIdx.y
#define LID2 threadIdx.z
#ifdef HIP
#define IMAGE3D hipTextureObject_t
#define IMAGE2D hipTextureObject_t
#else
#define IMAGE3D cudaTextureObject_t
#define IMAGE2D cudaTextureObject_t
#endif
#define MUINT2(a, b) make_uint2(a, b)
#ifndef LTYPE3
#define MINT3(a, b, c) make_int3(a, b, c)
#else
#define MINT3(a, b, c) make_long3(a, b, c)
#endif
#define MUINT3(a, b, c) make_uint3(a, b, c)
#define MFLOAT3(a, b, c) make_float3(a, b, c)
#define MFLOAT2(a, b) make_float2(a, b)
#define CMFLOAT3 make_float3
#define CMINT3 make_int3
#define CMINT4 make_int4
#define KERNEL extern "C" __global__
#define KERNEL2 KERNEL
#define KERNEL3 KERNEL
#define KERN KERNEL
#define BARRIER __syncthreads();
#define POWR __powf
#define POWN __powf
#define RCP(x) __frcp_rn(x)
#define FLOOR floorf
#define CEIL ceilf
#define ATAN2 atan2f
#define ACOS acosf
#define SQRT sqrtf
#define RSQRT(x) __frsqrt_rn(x)
#define LOG logf
#define CROSS cross
#define DISTANCE distance

template <typename tyT> 
inline __device__ tyT sign(tyT val) {
    return (tyT(0) < val) - (val < tyT(0));
}

inline __device__ float3 expf3(float3 a) {
	return make_float3(__expf(a.x), __expf(a.y), __expf(a.z));
}

inline __device__ float clamp(float f, float a, float b) {
    return fmaxf(a, fminf(f, b));
}

inline __device__ float3 clamp3(float3 a, float b, float3 c) {
	return make_float3(clamp(a.x, b, c.x), clamp(a.y, b, c.y), clamp(a.z, b, c.z));
}

inline __device__ float3 fdividef3(float3 a, float3 b) {
	return make_float3(fdividef(a.x, b.x), fdividef(a.y, b.y), fdividef(a.z, b.z));
}

inline __device__ int3 __float2int_rz3(float3 a) {
	return make_int3(__float2int_rz(a.x), __float2int_rz(a.y), __float2int_rz(a.z));
}

inline __device__ float3 make_int3_float3(int3 a) {
	return make_float3((float)a.x, (float)a.y, (float)a.z);
}

template <typename tyT> 
inline __device__ int3 make_long3_int3(tyT a) {
    return make_int3(static_cast<int>(a.x), static_cast<int>(a.y), static_cast<int>(a.z));
}

// CUDA's vector_types.h defines no arithmetic operators for int3/float2/float3, so we provide them
// here. HIP's HIP_vector_type already overloads all of these (+, -, *, /, unary -, compound
// assignment, and scalar/vector mixes), so defining them again makes every use ambiguous. Skip them
// under HIP and rely on the built-in ones. NOTE: the only behavioral difference is uint3 - int,
// which yields int3 here but uint3 under HIP; its single use is CFLOAT3(N - 1) where N >= 1, so the
// converted float value is identical.
#ifndef HIP
inline __device__ int3 operator-(const int3 a, const int3 b) {
	return make_int3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline __device__ int3 operator-(uint3 a, int b) {
	return make_int3(a.x - b, a.y - b, a.z - b);
}

inline __device__ int3 operator/(int3 a, int b) {
	return make_int3(a.x / b, a.y / b, a.z / b);
}

inline __device__ float2 operator-(float2 a, float2 b) {
	return make_float2(a.x - b.x, a.y - b.y);
}

inline __device__ int3 operator+(const int3 a, const int3 b) {
	return make_int3(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline __device__ float2 operator-(float a, float2 b) {
	return make_float2(a - b.x, a - b.y);
}

inline __device__ float2 operator+(float2 a, float2 b) {
	return make_float2(a.x + b.x, a.y + b.y);
}

inline __device__ float2 operator+(float a, float2 b) {
	return make_float2(a + b.x, a + b.y);
}

inline __device__ float2 operator+(float2 a, float b) {
	return make_float2(a.x + b, a.y + b);
}

inline __device__ float2 operator-(float2 a, float b) {
	return make_float2(a.x - b, a.y - b);
}

inline __device__ float2 operator*(float2 a, float2 b) {
	return make_float2(a.x * b.x, a.y * b.y);
}

inline __device__ float2 operator*(float b, float2 a) {
	return make_float2(a.x * b, a.y * b);
}

inline __device__ float2 operator*(float2 a, float b) {
	return make_float2(a.x * b, a.y * b);
}

inline __device__ float2 operator/(float2 a, float2 b) {
	return make_float2(a.x / b.x, a.y / b.y);
}

inline __device__ float2 operator/(float2 a, float b) {
	return make_float2(a.x / b, a.y / b);
}

inline __device__ float3 operator-(float3 a) {
	return make_float3(-a.x, -a.y, -a.z);
}

inline __device__ float3 operator-(float3 a, float3 b) {
	return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline __device__ float3 operator-(float3 a, float b) {
	return make_float3(a.x - b, a.y - b, a.z - b);
}

inline __device__ float3 operator+(float3 a, float3 b) {
	return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline __device__ float3 operator/(float3 a, float3 b) {
	return make_float3(a.x / b.x, a.y / b.y, a.z / b.z);
}

inline __device__ float3 operator/(float3 a, float b) {
	return make_float3(a.x / b, a.y / b, a.z / b);
}

inline __device__ float3 operator*(float3 a, float3 b) {
	return make_float3(a.x * b.x, a.y * b.y, a.z * b.z);
}

inline __device__ float3 operator*(float3 a, float b) {
	return make_float3(a.x * b, a.y * b, a.z * b);
}

inline __device__ float3 operator*(float b, float3 a) {
	return make_float3(a.x * b, a.y * b, a.z * b);
}

inline __device__ void operator*=(float3& a, float3 b) {
	a.x *= b.x;
	a.y *= b.y;
	a.z *= b.z;
}

inline __device__ void operator-=(float3& a, float3 b) {
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
}

inline __device__ void operator+=(float3& a, float3 b) {
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
}

inline __device__ void operator-=(float2& a, float b) {
	a.x -= b;
	a.y -= b;
}

inline __device__ void operator*=(float2& a, float b) {
	a.x *= b;
	a.y *= b;
}

inline __device__ void operator+=(float2& a, float2 b) {
	a.x += b.x;
	a.y += b.y;
}

inline __device__ void operator+=(float2& a, float b) {
	a.x += b;
	a.y += b;
}

inline __device__ void operator-=(float2& a, float2 b) {
	a.x -= b.x;
	a.y -= b.y;
}

inline __device__ void operator*=(float2& a, float2 b) {
	a.x *= b.x;
	a.y *= b.y;
}

inline __device__ void operator/=(float2& a, float b) {
	a.x /= b;
	a.y /= b;
}

inline __device__ void operator/=(float2& a, float2 b) {
	a.x /= b.x;
	a.y /= b.y;
}

inline __device__ void operator*=(float3& a, float b) {
	a.x *= b;
	a.y *= b;
	a.z *= b;
}

inline __device__ void operator-=(float3& a, float b) {
	a.x -= b;
	a.y -= b;
	a.z -= b;
}

inline __device__ void operator+=(float3& a, float b) {
	a.x += b;
	a.y += b;
	a.z += b;
}

inline __device__ void operator/=(float3& a, float b) {
	a.x /= b;
	a.y /= b;
	a.z /= b;
}

inline __device__ void operator/=(float3& a, float3 b) {
	a.x /= b.x;
	a.y /= b.y;
	a.z /= b.z;
}
#endif // !HIP (vector arithmetic operators are built into HIP_vector_type)

inline __device__ float3 fmin(float3 a, float3 b) {
    return make_float3(a.x < b.x ? a.x : b.x, a.y < b.y ? a.y : b.y, a.z < b.z ? a.z : b.z);
}

inline __device__ float3 fmax(float3 a, float3 b) {
    return make_float3(a.x > b.x ? a.x : b.x, a.y > b.y ? a.y : b.y, a.z > b.z ? a.z : b.z);
}

inline __device__ float dot(float3 a, float3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline __device__ float length(float3 v) {
    return sqrtf(dot(v, v));
}

inline __device__ float3 cross(float3 a, float3 b) {
    return make_float3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}

inline __device__ float3 normalize(float3 v) {
    float invLen = rsqrtf(dot(v, v));
    return v * invLen;
}

inline __device__ float3 __fmaf_rn3(float a, float3 b, float3 c) {
	return make_float3(__fmaf_rn(a, b.x, c.x), __fmaf_rn(a, b.y, c.y), __fmaf_rn(a, b.z, c.z));
}

inline __device__ float3 __fmaf_rn3(float3 a, float b, float3 c) {
	return make_float3(__fmaf_rn(a.x, b, c.x), __fmaf_rn(a.y, b, c.y), __fmaf_rn(a.z, b, c.z));
}

inline __device__ float3 __fmaf_rn3(float3 a, float3 b, float3 c) {
	return make_float3(__fmaf_rn(a.x, b.x, c.x), __fmaf_rn(a.y, b.y, c.y), __fmaf_rn(a.z, b.z, c.z));
}

inline __device__ float2 __fmaf_rn2(float a, float2 b, float2 c) {
	return make_float2(__fmaf_rn(a, b.x, c.x), __fmaf_rn(a, b.y, c.y));
}

inline __device__ float distance(float3 a, float3 b) {
    return length(a - b);
}

inline __device__ float3 fabs(float3 v) {
    return make_float3(fabs(v.x), fabs(v.y), fabs(v.z));
}

inline __device__ float2 fabs(float2 v) {
    return make_float2(fabs(v.x), fabs(v.y));
}

template<typename T>
inline __device__ float3 convert_float3(const T& a) {
    return make_float3(static_cast<float>(a.x), static_cast<float>(a.y), static_cast<float>(a.z));
}

template<typename T>
inline __device__ uint3 convert_uint3(const T& a) {
    return make_uint3(static_cast<unsigned int>(a.x), static_cast<unsigned int>(a.y), static_cast<unsigned int>(a.z));
}
#endif

#if STYPE == 1 || STYPE == 2 || STYPE == 4 || STYPE == 5
DEVICE void getIndex(int3* i, const uint d_size_x, const uint d_sizey, const uint currentSubset) {
#if STYPE == 1
	(*i).x *= NSUBSETS;
	(*i).x += currentSubset;
	(*i).y = (*i).x % d_sizey;
	(*i).z = (*i).x / (d_size_x * d_sizey);
	(*i).x /= d_sizey;
	(*i).x = (*i).x % d_size_x;
#elif STYPE == 2
	(*i).x *= NSUBSETS;
	(*i).x += currentSubset;
	(*i).z = (*i).x / (d_size_x * d_sizey);
	(*i).y = (*i).x / d_size_x;
	(*i).y = (*i).y % d_sizey;
	(*i).x = (*i).x % d_size_x;
#elif STYPE == 4
	(*i).x = ((*i).x / d_size_x) * d_size_x * NSUBSETS + d_size_x * currentSubset + (*i).x % d_size_x;
	(*i).z = (*i).x / (d_size_x * d_sizey);
	(*i).y = ((*i).x / d_size_x) % d_sizey;
	(*i).x = (*i).x % d_size_x;
#elif STYPE == 5
	(*i).y = (*i).x % d_sizey;
	(*i).x = ((*i).x / d_sizey * NSUBSETS + currentSubset);
	(*i).z = (*i).x / d_size_x;
	(*i).x = (*i).x % d_size_x;
#endif
}
#endif
#ifdef USEIMAGES
    #define IMTYPE IMAGE3D
    #ifdef MASKBP3D
        #define MASKBPTYPE IMAGE3D
    #else
        #define MASKBPTYPE IMAGE2D
    #endif
    #ifdef MASKFP3D
        #define MASKFPTYPE IMAGE3D
    #else
        #define MASKFPTYPE IMAGE2D
    #endif
#else
    #define IMTYPE const CLGLOBAL float* CLRESTRICT 
    #define MASKBPTYPE const CLGLOBAL uchar* CLRESTRICT 
    #define MASKFPTYPE const CLGLOBAL uchar* CLRESTRICT 

#endif

#if (defined(MASKBP) || defined(MASKBP3D) || defined(MASKPRIOR))
// Read the backprojection/prior mask value at the input voxel
// ind contains the voxel coordinates; for the ray-based projector type 4 these are normalized [0, 1) coordinates
// d_N contains the mask dimensions, only used with buffers (2D masks ignore ind.z)
DEVICE int readMaskBP(MASKBPTYPE maskBP,
	const T3 ind,
	const uint3 d_N) {
#ifdef USEIMAGES
#if defined(METAL)
#ifdef MASKBPNORM
    #ifdef MASKBP3D
        return static_cast<int>(metal::round(maskBP.sample(samplerMask, ind).r));
    #else
        return static_cast<int>(metal::round(maskBP.sample(samplerMask, ind.xy).r));
    #endif
#else
    #ifdef MASKBP3D
        return static_cast<int>(metal::round(maskBP.read(uint3(ind)).r));
    #else
        return static_cast<int>(metal::round(maskBP.read(uint2(ind.xy)).r));
    #endif
#endif
#elif defined(CUDA) || defined(HIP)
	return
    #ifdef MASKBP3D
        tex3D<unsigned char>(maskBP, ind.x, ind.y, ind.z);
    #else
        tex2D<unsigned char>(maskBP, ind.x, ind.y);
    #endif
    #else
	return
    #ifdef MASKBP3D
        read_imageui(maskBP, sampler_MASK, (T4)(ind.x, ind.y, ind.z, 0)).w;
    #else
        read_imageui(maskBP, sampler_MASK, (T2)(ind.x, ind.y)).w;
    #endif
#endif
#else
#ifdef MASKBPNORM
	// Normalized coordinates
	const LONG indX = CLONG_rtz(ind.x * CFLOAT(d_N.x));
	const LONG indY = CLONG_rtz(ind.y * CFLOAT(d_N.y));
#ifdef MASKBP3D
	const LONG indZ = CLONG_rtz(ind.z * CFLOAT(d_N.z));
#endif
#else
	const LONG indX = CLONG_rtz(ind.x);
	const LONG indY = CLONG_rtz(ind.y);
#ifdef MASKBP3D
	const LONG indZ = CLONG_rtz(ind.z);
#endif
#endif
	return maskBP[indX + indY * CLONG_rtz(d_N.x)
#ifdef MASKBP3D
		+ indZ * CLONG_rtz(d_N.x) * CLONG_rtz(d_N.y)
#endif
	];
#endif
}
#endif

#if defined(MASKFP)
DEVICE int readMaskFP(MASKFPTYPE maskFP, MASKFPT ind) {
	return
#ifdef USEIMAGES
#ifdef METAL
#ifdef MASKFP3D
        static_cast<int>(metal::round(maskFP.read(uint3(ind.x, ind.y, ind.z)).r));
#else
        static_cast<int>(metal::round(maskFP.read(uint2(ind.x, ind.y)).r));
#endif
#elif defined(CUDA) || defined(HIP)
#ifdef MASKFP3D
        tex3D<unsigned char>(maskFP, ind.x, ind.y, ind.z);
#else
        tex2D<unsigned char>(maskFP, ind.x, ind.y);
#endif
#else
#ifdef MASKFP3D
        read_imageui(maskFP, sampler_MASKFP, (int4)(ind.x, ind.y, ind.z, 0)).w;
#else
        read_imageui(maskFP, sampler_MASKFP, (int2)(ind.x, ind.y)).w;
#endif
#endif
#else
        maskFP[ind];
#endif
}
#endif


// This function was taken from: https://streamhpc.com/blog/2016-02-09/atomic-operations-for-floats-in-opencl-improved/
// Computes the atomic_add for floats
// NOTE: Includes code for an OpenCL extension that enables float atomics, but this is currently only supported by Intel and POCL
#if defined(ATOMICF) && !defined(ATOMIC) && !defined(ATOMIC32) && defined(OPENCL)
#if defined(NVIDIA)
void atomicAdd(__global float* p, float val)
{
    float prev;
    asm volatile(
        "atom.global.add.f32 %0, [%1], %2;" 
        : "=f"(prev) 
        : "l"(p) , "f"(val) 
        : "memory" 
    );
}
#elif defined(AMD)
void atomicAdd(volatile __global float* p, float val) {
    __asm__ volatile (
		"global_atomic_add_f32 %0, %1, off\n"
		: 
		: "v"(p), "v"(val) 
		: "memory"
    );
}
#elif defined(INTEL)
#pragma OPENCL EXTENSION cl_ext_float_atomics : enable
#define atomicAdd(a,b) atomic_fetch_add((volatile atomic_float *)(a),(b)) 
#else
void atomicAdd(volatile CLGLOBAL float *addr, float val) {
	union {
		unsigned int u32;
		float        f32;
	} next, expected, current;
	current.f32 = *addr;
	do {
		expected.f32 = current.f32;
		next.f32 = expected.f32 + val;
		current.u32 = atomic_cmpxchg((volatile CLGLOBAL unsigned int *)addr, expected.u32, next.u32);
	} while (current.u32 != expected.u32);
}
#endif
#endif

#ifdef TOF //////////////// TOF ////////////////
#define _2PI 0.3989423f

DEVICE float normPDF(const float x, const float mu, const float invSigma, const float piPerSigma) {
	const float a = (x - mu) * invSigma;
	return piPerSigma * EXP(-0.5f * a * a);
}

DEVICE void TOFDis(const float3 diff, const float tc, const float LL, float* D, float* DD) {
	*D = length(diff * tc) - LL / 2.f;
	*DD = *D;
}

DEVICE float TOFWeight(const float element, const float invSigma, const float piPerSigma, const float D, const float DDsign, const float TOFCenter, float dX) {
	float output = normPDF(D, TOFCenter, invSigma, piPerSigma);
	dX *= DDsign;
#pragma unroll
	for (int tr = 1; tr < CINT(TRAPZ_BINS) - 1; tr++)
#ifdef USEMAD
		output += (normPDF(FMAD(-dX, CFLOAT(tr), D), TOFCenter, invSigma, piPerSigma) * 2.f);
	output += normPDF(FMAD(-element, DDsign, D), TOFCenter, invSigma, piPerSigma);
#else
		output += (normPDF(D - dX * CFLOAT(tr), TOFCenter, invSigma, piPerSigma) * 2.f);
	output += normPDF(D - element * DDsign, TOFCenter, invSigma, piPerSigma);
#endif
	return output;
}


DEVICE float TOFLoop(const float DDsign, const float element, CONSTANT float* TOFCenter, const float invSigma, const float piPerSigma, float* D, const float epps, float* TOFWeights) {
	float TOFSum = 0.f;
	const float dX = element / (TRAPZ_BINS - 1.f);
#if !defined(__CUDACC__) && !defined(__HIPCC__)
#pragma unroll NBINS
#endif
	for (int to = 0; to < NBINS; to++) {
		TOFWeights[to] = TOFWeight(element, invSigma, piPerSigma, *D, DDsign, TOFCenter[to], dX) * dX;
		TOFSum += TOFWeights[to];
	}
	if (TOFSum < epps)
		TOFSum = epps;
	return TOFSum;
}
#endif //////////////// END TOF ////////////////

#if defined(N_RAYS)
DEVICE void multirayCoordinateShiftXY(PTR_THR FLOAT3 *s, PTR_THR FLOAT3 *d, const int lor, const float cr) {
	float interval = cr / (CFLOAT(N_RAYS2D * 2));
	(*s).x += (interval - cr / 2.f);
	(*d).x += (interval - cr / 2.f);
	(*s).y += (interval - cr / 2.f);
	(*d).y += (interval - cr / 2.f);
	interval *= 2.f;
	(*s).x += interval * lor;
	(*d).x += interval * lor;
	(*s).y += interval * lor;
	(*d).y += interval * lor;
}

DEVICE void multirayCoordinateShiftZ(PTR_THR FLOAT3 *s, PTR_THR FLOAT3 *d, const int lor, const float cr) {
	float interval = cr / (CFLOAT(N_RAYS3D * 2));
	(*s).z += (interval - cr / 2.f);
	(*d).z += (interval - cr / 2.f);
	interval *= 2.f;
	(*s).z += interval * lor;
	(*d).z += interval * lor;
}
#endif

#if defined(FP) && !defined(PROJ5)
// Computes the forward projection
// Separate cases for the Siddon and interpolated projectors
DEVICE void forwardProject(const float local_ele, PTR_THR float *ax, const typeT local_ind, IMTYPE d_OSEM) {
#if defined(CUDA) || defined(HIP)
#ifdef USEIMAGES
    // if (local_ind.x <= 1.f && local_ind.y <= 1.f && local_ind.z <= 1.f && local_ind.x >= 0.f && local_ind.y >= 0.f && local_ind.z >= 0.f)
		*ax = (local_ele * tex3D<float>(d_OSEM, local_ind.x, local_ind.y, local_ind.z));
#else
	*ax = (local_ele * d_OSEM[local_ind]);
#endif
#elif defined(OPENCL)
#ifdef PTYPE4
    // if (local_ind.x <= 1.f && local_ind.y <= 1.f && local_ind.z <= 1.f && local_ind.x >= 0.f && local_ind.y >= 0.f && local_ind.z >= 0.f)
		*ax = (local_ele * read_imagef(d_OSEM, samplerForw, (T4)(local_ind, (typeTT)0)).w);
#else
#ifdef USEIMAGES
	*ax = (local_ele * read_imagef(d_OSEM, samplerSiddon, (T4)(local_ind, (typeTT)0)).w);
#else
	*ax = (local_ele * d_OSEM[local_ind]);
#endif
#endif
#elif defined(METAL)
#ifdef PTYPE4
    *ax = local_ele * d_OSEM.sample(samplerForw, (local_ind)).r;
#elif defined(USEIMAGES)
    *ax = local_ele * d_OSEM.read(uint3(local_ind.x, local_ind.y, local_ind.z)).r;
#else
    *ax = (local_ele * d_OSEM[local_ind]);
#endif
#endif
}

// Computes the forward projection
// Includes TOF-specific weighting
DEVICE void denominator(PTR_THR float *ax, const typeT localInd, float local_ele, IMTYPE d_OSEM
#ifdef TOF
	, const float TOFSum, float* TOFWeights
#ifdef LISTMODE
	, const int TOFIndex
#endif
#endif
) {
	float apu = 0.f;
	forwardProject(local_ele, &apu, localInd, d_OSEM);
#ifdef TOF
	const float dX = apu / TOFSum;
#if defined(LISTMODE) && !defined(SENS)
	int to = TOFIndex;
#else
#if !defined(__CUDACC__) && !defined(__HIPCC__)
#pragma unroll NBINS
#endif
	for (int to = 0; to < NBINS; to++) {
#endif
		const float joku = TOFWeights[to] * dX;
		ax[to] += joku;
#if !defined(LISTMODE) || defined(SENS)
	}
#endif
#else
	ax[0] += apu;
#endif
}
#endif

#if defined(BP) && !defined(PROJ5) && (defined(ATOMIC) || defined(ATOMIC32) || defined(ATOMICF))
// Compute the backprojection
DEVICE void rhs(const float local_ele, PTR_THR const float *ax, const LONG local_ind, CLGLOBAL CAST* d_rhs_OSEM, const uchar no_norm, CLGLOBAL CAST* d_Summ
#ifdef TOF
	, const FLOAT TOFSum, float* TOFWeights
#ifdef LISTMODE
	, const int TOFIndex
#endif
#endif
) {
#ifdef TOF
	FLOAT val = FLOAT_ZERO;
	const FLOAT dX = local_ele / TOFSum;
	FLOAT yaxTOF = FLOAT_ZERO;
#if defined(LISTMODE) && !defined(SENS)
	int to = TOFIndex;
#else
#if !defined(__CUDACC__) && !defined(__HIPCC__)
#pragma unroll NBINS
#endif
	for (int to = 0; to < NBINS; to++) {
#endif
		const float apu = dX * TOFWeights[to];
		val += apu;
		yaxTOF += (apu * ax[to]);
#if !defined(LISTMODE) || defined(SENS)
	}
#endif
#else
	float yaxTOF = ax[0] * local_ele;
	const float val = local_ele;
	// test
	if (ISNAN(local_ele))
		yaxTOF = 0.f;
	//else
	//	yaxTOF = 1.f;

#endif
#ifdef ATOMIC
	atom_add(&d_rhs_OSEM[local_ind], convert_long(yaxTOF * TH));
#elif defined(ATOMIC32)
#ifdef METAL
	atomicAdd((volatile device metal::atomic_int*)(&d_rhs_OSEM[local_ind]), CINT_rtz(yaxTOF * TH));
#else
	atomic_add(&d_rhs_OSEM[local_ind], CINT(yaxTOF * TH));
#endif
#else
#ifdef METAL
	atomicAdd((volatile device metal::atomic_float*)(&d_rhs_OSEM[local_ind]), yaxTOF);
#else
	atomicAdd((&d_rhs_OSEM[local_ind]), yaxTOF);
#endif
#endif
	if (no_norm == 0u)
#ifdef ATOMIC
		atom_add(&d_Summ[local_ind], convert_long(val * TH));
#elif defined(ATOMIC32)
#ifdef METAL
		atomicAdd((volatile device metal::atomic_int*)(&d_Summ[local_ind]), CINT_rtz(val * TH));
#else
		atomic_add(&d_Summ[local_ind], CINT(val * TH));
#endif
#else
#ifdef METAL
		atomicAdd((volatile device metal::atomic_float*)&d_Summ[local_ind], val);
#else
		atomicAdd((&d_Summ[local_ind]), val);
#endif
#endif
}
#endif


// Detector coordinates for listmode data
#ifdef LISTMODE
#ifdef INDEXBASED
DEVICE void getDetectorCoordinatesListmode(
#if defined(USEGLOBAL)
	const CLGLOBAL float* d_xy, const CLGLOBAL float* d_z, 
#else
	CONSTANT float* d_xy, CONSTANT float* d_z, 
#endif
	const CLGLOBAL ushort* trIndex, const CLGLOBAL ushort* axIndex, float3* s, float3* d, const size_t idx
#if defined(N_RAYS)
	, const int lorXY, const int lorZ, const float2 cr
#endif
) {
	const size_t i = idx * 2;
	size_t id = trIndex[i] * 2;
	size_t idz = axIndex[i];
	*s = CMFLOAT3(d_xy[id], d_xy[id + 1], d_z[idz]);
	id = trIndex[i + 1] * 2;
	idz = axIndex[i + 1];
	*d = CMFLOAT3(d_xy[id], d_xy[id + 1], d_z[idz]);
#if defined(N_RAYS)
	if (N_RAYS3D > 1)
		multirayCoordinateShiftZ(s, d, lorZ, cr.y);
	if (N_RAYS2D > 1)
		multirayCoordinateShiftXY(s, d, lorXY, cr.x);
#endif
}
#else
DEVICE void getDetectorCoordinatesListmode(const CLGLOBAL float* d_xyz, float3* s, float3* d, const size_t idx
#if defined(N_RAYS)
	, const int lorXY, const int lorZ, const float2 cr
#endif
) {
	const size_t i = idx * 6;
	*s = CMFLOAT3(d_xyz[i], d_xyz[i + 1], d_xyz[i + 2]);
	*d = CMFLOAT3(d_xyz[i + 3], d_xyz[i + 4], d_xyz[i + 5]);
#if defined(N_RAYS)
	if (N_RAYS3D > 1)
		multirayCoordinateShiftZ(s, d, lorZ, cr.y);
	if (N_RAYS2D > 1)
		multirayCoordinateShiftXY(s, d, lorXY, cr.x);
#endif
}
#endif
#endif

// Detector coordinates for CT data
#if defined(CT) && !defined(SPECTMASK)
#if !defined(USEGLOBAL)
DEVICE void getDetectorCoordinatesCT(CONSTANT float* d_xyz, 
#else
DEVICE void getDetectorCoordinatesCT(const CLGLOBAL float* CLRESTRICT d_xyz, 
#endif
#ifdef HELICAL
	const float r, 
#endif
#if !defined(USEGLOBAL)
	CONSTANT float* d_uv, 
#else
	const CLGLOBAL float* CLRESTRICT d_uv, 
#endif
	PTR_THR float3* s, PTR_THR float3* d, const int3 i, const uint d_size_x, const uint d_sizey, const float2 d_dPitch
#ifdef PROJ5
	, PTR_THR float3* dR, PTR_THR float3* dL, PTR_THR float3* dU, PTR_THR float3* dD
#endif
) {
	int id = i.z * 6;
	*s = CMFLOAT3(d_xyz[id], d_xyz[id + 1], d_xyz[id + 2]);
	*d = CMFLOAT3(d_xyz[id + 3], d_xyz[id + 4], d_xyz[id + 5]);
	const float2 indeksi = MFLOAT2(CFLOAT(i.x) - CFLOAT(d_size_x) / 2.f + .5f, CFLOAT(i.y) - CFLOAT(d_sizey) / 2.f + .5f);
#ifdef HELICAL
	const float angle = d_uv[i.z];
	const float dtheta = (d_dPitch.x / r) * indeksi.x;
	(*d).x += r * (COSF(angle + dtheta) - COSF(angle));
	(*d).y -= r * (SINF(angle + dtheta) - SINF(angle));
	(*d).z += indeksi.y * d_dPitch.y;
#else
	id = i.z * NA;
#if defined(PITCH)
	const float3 apuX = MFLOAT3(d_uv[id], d_uv[id + 1], d_uv[id + 2]);
	const float3 apuY = MFLOAT3(d_uv[id + 3], d_uv[id + 4], d_uv[id + 5]);
#ifdef USEMAD
	*d += FMAD3(apuX, indeksi.x, apuY * indeksi.y);
#else
	*d += apuX * indeksi.x + apuY * indeksi.y;
#endif
#ifdef PARALLEL
#ifdef USEMAD
	*s += FMAD3(apuX, indeksi.x, apuY * indeksi.y);
#else
	*s += apuX * indeksi.x + apuY * indeksi.y;
#endif
#endif
#if defined(PROJ5) && defined(FP)
#ifdef USEMAD
	*dR = FMAD3(-apuX, 0.5f, *d);
	*dL = FMAD3(apuX, 0.5f, *d);
	*dU = FMAD3(apuY, 0.5f, *d);
	*dD = FMAD3(-apuY, 0.5f, *d);
#else
	*dR = *d - apuX * 0.5f;
	*dL = *d + apuX * 0.5f;
	*dU = *d + apuY * 0.5f;
	*dD = *d - apuY * 0.5f;
#endif
#endif
#else
	const float apuX = d_uv[id];
	const float apuY = d_uv[id + 1];
	(*d).x += indeksi.x * apuX;
	(*d).y += indeksi.x * apuY;
	(*d).z += indeksi.y * d_dPitch.y;
#ifdef PARALLEL
	(*s).x += indeksi.x * apuX;
	(*s).y += indeksi.x * apuY;
	(*s).z += indeksi.y * d_dPitch.y;
#endif
#if defined(PROJ5) && defined(FP)
#ifdef USEMAD
	*dR = CMFLOAT3(FMAD(-apuX, 0.5f, (*d).x), FMAD(-apuY, 0.5f, (*d).y), (*d).z);
	*dL = CMFLOAT3(FMAD(apuX, 0.5f, (*d).x), FMAD(apuY, 0.5f, (*d).y), (*d).z);
	*dU = CMFLOAT3((*d).x, (*d).y, FMAD(d_dPitch.y, 0.5f, (*d).z));
	*dD = CMFLOAT3((*d).x, (*d).y, FMAD(-d_dPitch.y, 0.5f, (*d).z));
#else
	*dR = CMFLOAT3((*d).x - apuX * 0.5f, (*d).y - apuY * 0.5f, (*d).z);
	*dL = CMFLOAT3((*d).x + apuX * 0.5f, (*d).y + apuY * 0.5f, (*d).z);
	*dU = CMFLOAT3((*d).x, (*d).y, (*d).z + d_dPitch.y * 0.5f);
	*dD = CMFLOAT3((*d).x, (*d).y, (*d).z - d_dPitch.y * 0.5f);
#endif
#endif
#endif
#endif
}

#elif defined(SPECT)
DEVICE void extendRayToEllipse(
    PTR_THR FLOAT3 *s, // Ray start point
    PTR_THR FLOAT3 *d, // Ray end point
    const FLOAT3 ellipseCenter,
    const FLOAT3 ellipseRadii,
    const FLOAT ellipsePower
) {
#ifdef CUPY_HIP_FINITE_ELLIPSE_POWER
    // CuPy stages the public infinity sentinel as a finite value; this
    // comparison is safe under hipRTC's finite-math assumptions.
    if (ellipsePower > 1.0e20f) {
#else
    if (ISINF(ellipsePower)) {
#endif
        const FLOAT cx = ellipseCenter.x;
        const FLOAT cy = ellipseCenter.y;
        const FLOAT cz = ellipseCenter.z;
        const FLOAT ax = ellipseRadii.x;
        const FLOAT ay = ellipseRadii.y;
        const FLOAT az = ellipseRadii.z;
        const FLOAT3 boxMin = CMFLOAT3(cx - ax, cy - ay, cz - az);
        const FLOAT3 boxMax = CMFLOAT3(cx + ax, cy + ay, cz + az);
        const FLOAT3 p0 = *s;
        const FLOAT3 p1 = *d;
        const FLOAT3 dir = p1 - p0;

        FLOAT tmin = -1e8f;
        FLOAT tmax =  1e8f;
        const FLOAT epsVal = 1.0e-8f;

        if (FABS(dir.x) < epsVal) {
            if (p0.x < boxMin.x || p0.x > boxMax.x) {
                *d = *s;
                return;
            }
        } else {
            FLOAT t1 = (boxMin.x - p0.x) / dir.x;
            FLOAT t2 = (boxMax.x - p0.x) / dir.x;
            FLOAT tNear = FMIN(t1, t2);
            FLOAT tFar = FMAX(t1, t2);
            tmin = FMAX(tmin, tNear);
            tmax = FMIN(tmax, tFar);
        }

        if (FABS(dir.y) < epsVal) {
            if (p0.y < boxMin.y || p0.y > boxMax.y) {
                *d = *s;
                return;
            }
        } else {
            FLOAT t1 = (boxMin.y - p0.y) / dir.y;
            FLOAT t2 = (boxMax.y - p0.y) / dir.y;
            FLOAT tNear = FMIN(t1, t2);
            FLOAT tFar = FMAX(t1, t2);
            tmin = FMAX(tmin, tNear);
            tmax = FMIN(tmax, tFar);
        }

        if (FABS(dir.z) < epsVal) {
            if (p0.z < boxMin.z || p0.z > boxMax.z) {
                *d = *s;
                return;
            }
        } else {
            FLOAT t1 = (boxMin.z - p0.z) / dir.z;
            FLOAT t2 = (boxMax.z - p0.z) / dir.z;
            FLOAT tNear = FMIN(t1, t2);
            FLOAT tFar = FMAX(t1, t2);
            tmin = FMAX(tmin, tNear);
            tmax = FMIN(tmax, tFar);
        }

        if (tmax < tmin) {
            *d = *s;
            return;
        }

        if (!((p0.x >= boxMin.x && p0.x <= boxMax.x) && (p0.y >= boxMin.y && p0.y <= boxMax.y) && (p0.z >= boxMin.z && p0.z <= boxMax.z)))
            *s = p0 + tmin * dir;

        *d = p0 + tmax * dir;
    } else if (ellipsePower == FLOAT_TWO) {
        const FLOAT cx = ellipseCenter.x;
        const FLOAT cy = ellipseCenter.y;
        const FLOAT cz = ellipseCenter.z;
        const FLOAT ax = ellipseRadii.x;
        const FLOAT ay = ellipseRadii.y;
        const FLOAT az = ellipseRadii.z;
        const FLOAT3 p0 = *s;
        const FLOAT3 dir = *d - p0;
        const FLOAT px = p0.x - cx;
        const FLOAT py = p0.y - cy;
        const FLOAT invAx2 = FLOAT_ONE / (ax * ax);
        const FLOAT invAy2 = FLOAT_ONE / (ay * ay);
        const FLOAT A = dir.x * dir.x * invAx2 + dir.y * dir.y * invAy2;
        const FLOAT B = FLOAT_TWO * (px * dir.x * invAx2 + py * dir.y * invAy2);
        const FLOAT C = px * px * invAx2 + py * py * invAy2 - FLOAT_ONE;
        const FLOAT epsVal = 1.0e-8f;
        FLOAT tmin = -1e8f;
        FLOAT tmax = 1e8f;

        if (A < epsVal) {
            if (C > FLOAT_ZERO) {
                *d = *s;
                return;
            }
        } else {
            const FLOAT discriminant = B * B - 4.f * A * C;
            if (discriminant < FLOAT_ZERO) {
                *d = *s;
                return;
            }
            const FLOAT t1 = (-B - SQRT(discriminant)) / (FLOAT_TWO * A);
            const FLOAT t2 = (-B + SQRT(discriminant)) / (FLOAT_TWO * A);
            tmin = FMAX(tmin, t1);
            tmax = FMIN(tmax, t2);
        }

        if (FABS(dir.z) < epsVal) {
            if (p0.z < cz - az || p0.z > cz + az) {
                *d = *s;
                return;
            }
        } else {
            const FLOAT t1 = (cz - az - p0.z) / dir.z;
            const FLOAT t2 = (cz + az - p0.z) / dir.z;
            tmin = FMAX(tmin, FMIN(t1, t2));
            tmax = FMIN(tmax, FMAX(t1, t2));
        }

        if (tmax < tmin) {
            *d = *s;
            return;
        }

        if (C > FLOAT_ZERO || p0.z < cz - az || p0.z > cz + az)
            *s = p0 + tmin * dir;

        *d = p0 + tmax * dir;
    }
}

// SPECT sinogram coordinates
DEVICE void getDetectorCoordinatesSPECT(
#if defined(USEGLOBAL)
	const CLGLOBAL float *d_xyz,
    const CLGLOBAL float *d_uv, 
#else
	CONSTANT float *d_xyz,
    CONSTANT float *d_uv, 
#endif
    PTR_THR FLOAT3 *s, // Ray start point
    PTR_THR FLOAT3 *d, // Ray end point
    const int3 i,
    const uint d_size_x, // Detector element count x-direction
    const uint d_sizey, // Detector element count y-direction
    const FLOAT2 d_dPitch, // Detector element size [mm]
    const CLGLOBAL float* d_rayShiftsDetector, // Ray shifts [mm]
    const CLGLOBAL float* d_rayShiftsSource, // Ray shifts [mm]
    const CLGLOBAL uint* d_detectorVector, // Detector head for each projection
    int lor,
    const FLOAT3 ellipseCenter,
    const FLOAT3 ellipseRadii,
    const FLOAT ellipsePower
#if defined(ORTH)
    , PTR_THR FLOAT3 *collimatorOrigin
#endif
) {
	uint id = i.z * 6;
	*s = CMFLOAT3((FLOAT)d_xyz[id], (FLOAT)d_xyz[id + 1], (FLOAT)d_xyz[id + 2]); // TODO remove cast
	*d = CMFLOAT3((FLOAT)d_xyz[id + 3], (FLOAT)d_xyz[id + 4], (FLOAT)d_xyz[id + 5]); // TODO remove cast
	const FLOAT2 shift_det_elem = MFLOAT2(
        d_dPitch.x * (CFLOAT(i.x) + (FLOAT_ONE - CFLOAT(d_size_x)) * FLOAT_HALF),
        d_dPitch.y * (CFLOAT(i.y) + (FLOAT_ONE - CFLOAT(d_sizey)) * FLOAT_HALF)
    ); // Amount of shift from sinogram center to current detector element
	
    id = i.z * NA; // Index of d_uv (detector panel normal vector)
    const uint detectorElement = i.x + i.y * d_size_x;
    const uint detectorHead = d_detectorVector[i.z];
    uint idShift = 2*lor + (2*N_RAYS) * (detectorElement + detectorHead * d_size_x * d_sizey); // Index of rayShiftsDetector

	const FLOAT apuX = d_uv[id]; // X component of detector panel normal vector
	const FLOAT apuY = d_uv[id + 1]; // Y component of detector panel normal vector
    
	(*d).x += apuX * (shift_det_elem.x + d_rayShiftsDetector[idShift]); // Shift to current element + shift to rayShiftsDetector
	(*d).y += apuY * (shift_det_elem.x + d_rayShiftsDetector[idShift]);
	(*d).z += shift_det_elem.y + d_rayShiftsDetector[idShift+1];
	(*s).x += apuX * (shift_det_elem.x + d_rayShiftsSource[idShift]);
	(*s).y += apuY * (shift_det_elem.x + d_rayShiftsSource[idShift]);
	(*s).z += shift_det_elem.y + d_rayShiftsSource[idShift+1];

    #if defined(ORTH)
    *collimatorOrigin = *s;
    #endif
    extendRayToEllipse(s, d, ellipseCenter, ellipseRadii, ellipsePower);
}
#else
#if defined(RAW) || defined(SENS)
// Get the detector coordinates for the current (raw) measurement
DEVICE void getDetectorCoordinatesRaw(
#if defined(USEGLOBAL)
	const CLGLOBAL float* d_xy, const CLGLOBAL float* d_z,
#else
	CONSTANT float* d_xy, CONSTANT float* d_z,
#endif
	const int3 i, float3* s, float3* d, const int2 indz
#if defined(N_RAYS)
	, const int lorXY, const int lorZ, const float2 cr
#endif
) {
	(*s).x = d_xy[i.x * 2];
	(*s).y = d_xy[i.x * 2 + 1];
	(*d).x = d_xy[i.y * 2];
	(*d).y = d_xy[i.y * 2 + 1];
	(*s).z = d_z[indz.x];
	(*d).z = d_z[indz.y];
#if defined(N_RAYS)
	if (N_RAYS3D > 1)
		multirayCoordinateShiftZ(s, d, lorZ, cr.y);
	if (N_RAYS2D > 1)
		multirayCoordinateShiftXY(s, d, lorXY, cr.x);
#endif
}
#endif


#if !defined(RAW) && !defined(LISTMODE) && !defined(CT) && !defined(SPECT) && !defined(PET)
// Get the detector coordinates for the current sinogram bin (index-based subsets)
DEVICE void getDetectorCoordinates(const CLGLOBAL uint *d_xyindex, const CLGLOBAL ushort *d_zindex, const size_t idx,
	PTR_THR float3* s, PTR_THR float3* d, 
#if !defined(USEGLOBAL)
	CONSTANT float *d_xy, CONSTANT float *d_z
#else
	const CLGLOBAL float *d_xy, const CLGLOBAL float *d_z
#endif
#if defined(N_RAYS)
	, const int lorXY, const int lorZ, const float2 cr
#endif
#if defined(NLAYERS)
	, const uint d_sizey, const uint d_size_x
#endif
) {
	const uint ind = d_xyindex[idx] * 4;
#if defined(NLAYERS)
	const uint indz = d_zindex[idx] * 3;
	const int layer = CINT(d_z[indz]);
	(*s).x = d_xy[ind + layer * d_size_x * d_sizey];
	(*s).y = d_xy[ind + 1 + layer * d_size_x * d_sizey];
	(*d).x = d_xy[ind + 2 + layer * d_size_x * d_sizey];
	(*d).y = d_xy[ind + 3 + layer * d_size_x * d_sizey];
	(*s).z = d_z[indz + 1];
	(*d).z = d_z[indz + 2];
#else
	const uint indz = d_zindex[idx] * 2;
	(*s).x = d_xy[ind];
	(*s).y = d_xy[ind + 1];
	(*d).x = d_xy[ind + 2];
	(*d).y = d_xy[ind + 3];
	(*s).z = d_z[indz];
	(*d).z = d_z[indz + 1];
#endif
#if defined(N_RAYS)
	if (N_RAYS3D > 1)
		multirayCoordinateShiftZ(s, d, lorZ, cr.y);
	if (N_RAYS2D > 1)
		multirayCoordinateShiftXY(s, d, lorXY, cr.x);
#endif
}
#endif
#endif

#if !defined(SUBSETS) && !defined(CT)
// Get the detector coordinates for the current measurement (no subsets or using full sinogram subsets)
DEVICE void getDetectorCoordinatesFullSinogram(const uint d_size_x, const int3 i, PTR_THR FLOAT3* s, PTR_THR FLOAT3* d, 
#if defined(USEGLOBAL)
	const CLGLOBAL float* d_xy, const CLGLOBAL float* d_z
#else
	CONSTANT float* d_xy, CONSTANT float* d_z
#endif
#if defined(N_RAYS)
	, const int lorXY, const int lorZ, const float2 cr
#endif
#if defined(NLAYERS)
	, const uint d_sizey, const uint layer
#endif
) {
	const int id = (i.x + i.y * d_size_x) * 4;
	const int idz = i.z * 2;
#if defined(NLAYERS)
	*s = CMFLOAT3(d_xy[id + layer * d_size_x * d_sizey * 4], d_xy[id + layer * d_size_x * d_sizey * 4 + 1], d_z[idz]);
	*d = CMFLOAT3(d_xy[id + layer * d_size_x * d_sizey * 4 + 2], d_xy[id + layer * d_size_x * d_sizey * 4 + 3], d_z[idz + 1]);
#else
	*s = CMFLOAT3(d_xy[id], d_xy[id + 1], d_z[idz]);
	*d = CMFLOAT3(d_xy[id + 2], d_xy[id + 3], d_z[idz + 1]);
#endif
#if defined(N_RAYS)
	if (N_RAYS3D > 1)
		multirayCoordinateShiftZ(s, d, lorZ, cr.y);
	if (N_RAYS2D > 1)
		multirayCoordinateShiftXY(s, d, lorXY, cr.x);
#endif
}
#endif

#if defined(ATN) && !defined(CT)
DEVICE void compute_attenuation(const float val, const typeT ind, IMTYPE d_atten, PTR_THR float *jelppi, const int ii) {
	if (ii == 0) {
        *jelppi += val * -
#if defined(CUDA) || defined(HIP)
#ifdef USEIMAGES
		tex3D<float>(d_atten, ind.x, ind.y, ind.z);
#else
		d_atten[ind];
#endif
#elif defined(METAL)
#if defined(PTYPE4)
        d_atten.sample(samplerForw, float3(ind.x, ind.y, ind.z)).r;
#elif defined(USEIMAGES)
        d_atten.read(uint3(ind.x, ind.y, ind.z)).r;
#else
        d_atten[ind];
#endif
#else
#if defined(PTYPE4)
		read_imagef(d_atten, samplerForw, (float4)(ind.x, ind.y, ind.z, 0.f)).w;
#else
#ifdef USEIMAGES
		read_imagef(d_atten, samplerSiddon, (int4)(ind.x, ind.y, ind.z, 0)).w;
#else
		d_atten[ind];
#endif
#endif
#endif
	}
}
#endif

#if !defined(PTYPE4) && !defined(PROJ5)
// Compute the voxel index where the current perpendicular measurement starts
DEVICE int perpendicular_start(const float d_b, const float d, const float d_d, const uint d_N) {
	int tempi = CINT(FLOOR((d - d_b) / d_d));
	if (tempi < 0 || tempi >= CINT(d_N))
		tempi = 0;
	return tempi;
}

// Compute the probability for the perpendicular elements
DEVICE void perpendicular_elements(const float d_b, const float d_d1, const uint d_N1, const float d, const float d_d2, const uint d_N2, 
	PTR_THR float* templ_ijk, PTR_THR int3* tempi, PTR_THR LONG *z_loop, const uint d_N, const uint d_NN, 
	const size_t idx, const FLOAT global_factor, const FLOAT local_scat, 
#if !defined(CT) && defined(ATN) && !defined(ATNM)
	IMTYPE d_atten, const int ii, 
#elif !defined(CT) && !defined(ATN) && defined(ATNM)
	const CLGLOBAL float* CLRESTRICT d_atten,
#endif
	const float local_norm, const float L) {
	int apu = perpendicular_start(d_b, d, d_d1, d_N1);
	*z_loop = CLONG_rtz(apu) * CLONG_rtz(d_N) + *z_loop * CLONG_rtz(d_N1) * CLONG_rtz(d_N2);
	if (d_N == 1)
		(*tempi).x = apu;
	else
		(*tempi).y = apu;
#ifdef CT //////////////// CT ////////////////
#ifdef N_RAYS //////////////// MULTIRAY ////////////////
	* templ_ijk = 1.f / CFLOAT(N_RAYS);
#else
	* templ_ijk = 1.f;
#endif //////////////// END MULTIRAY ////////////////
#else //////////////// PET ////////////////
	// Probability
#ifdef N_RAYS //////////////// MULTIRAY ////////////////
#ifdef TOTLENGTH
	float temp = FLOAT_ONE / (L * CFLOAT(N_RAYS));
#else
	float temp = FLOAT_ONE / (CFLOAT(d_N2) * d_d2 * CFLOAT(N_RAYS));
#endif
#elif defined(ORTH)
	float temp = FLOAT_ONE;
#else
#ifdef TOTLENGTH
	float temp = FLOAT_ONE / L;
#else
	float temp = FLOAT_ONE / (CFLOAT(d_N2) * d_d2);
#endif
#endif //////////////// END MULTIRAY ////////////////
#if defined(ATN) && !defined(SPECT) //////////////// ATTENUATION ////////////////
		float jelppi = FLOAT_ZERO;
#ifdef USEIMAGES
		int3 atnind = *tempi;
#else
		LONG atnind = *z_loop;
#endif
		for (int iii = 0u; iii < d_N2; iii++) {
#ifdef USEIMAGES
			if (d_NN == 1)
				atnind.x = iii;
			else
				atnind.y = iii;
#else
			// The image branch above sets the varying coordinate, the linear index has to do the same
		// (*z_loop has the varying coordinate at zero) instead of accumulating the offset
		atnind = *z_loop + CLONG_rtz(iii) * CLONG_rtz(d_NN);
#endif
			compute_attenuation(d_d2, atnind, d_atten, &jelppi, ii);
		}
		temp *= EXP(jelppi);
#endif //////////////// END ATTENUATION ////////////////
#ifdef NORM
		temp *= local_norm;
#endif
#ifdef SCATTER
		temp *= local_scat;
#endif
#ifdef ATNM
		temp *= d_atten[idx];
#endif
	temp *= global_factor;
	*templ_ijk = temp;
#endif //////////////// END PET OR CT ////////////////
}
#endif

#if defined(SIDDON)
// Compute functions (9) and (29) (detector larger than source)
DEVICE void d_g_s_precomp(const float tmin, const float t_min, const float tmax, const float t_max, PTR_THR uint *v_min, PTR_THR uint *v_max, PTR_THR float *t_0, PTR_THR int *v_u, 
	const float diff, const float b, const float d, const float s, const uint N) {

	if (tmin == t_min)
		// (11)
		*v_min = 1u;
	else {
		// (2) and (19)
		const float p_t = s + tmin * (diff);
		// (12)
		*v_min = CUINT_rtp((p_t - b) / d);
	}
	if (tmax == t_max)
		// (13)
		*v_max = N;
	else {
		// (2) and (19)
		const float p_t = s + tmax * (diff);
		// (14)
		*v_max = CUINT_sat_rtz((p_t - b) / d);
	}
	// (9)
	*t_0 += ((CFLOAT(*v_min) * d) / (diff));
	//  (29)
	*v_u = 1;
}

// Compute functions (9) and (29) (source larger than detector)
DEVICE void s_g_d_precomp(const float tmin, const float t_min, const float tmax, const float t_max, PTR_THR uint *v_min, PTR_THR uint *v_max, PTR_THR float *t_0, PTR_THR int *v_u, 
	const float diff, const float b, const float d, const float s, const uint N) {

	if (tmin == t_min)
		// (15)
		*v_max = N - 1u;
	else {
		// (2) and (19)
		const float p_t = s + tmin * (diff);
		// (16)
		*v_max = CUINT_sat_rtz((p_t - b) / d);
	}
	if (tmax == t_max)
		// (17)
		*v_min = 0u;
	else {
		// (2) and (19)
		const float p_t = s + tmax * (diff);
		// (18)
		*v_min = CUINT_rtp((p_t - b) / d);
	}
	// (9)
	*t_0 += ((CFLOAT(*v_max) * d) / (diff));
	// (29)
	*v_u = -1;
}

// Compute the index of the current voxel
DEVICE LONG compute_ind(const int tempj, const int tempi, const int tempk, const uint d_Nx, const uint d_Nyx) {
	LONG local_ind = CLONG_rtz(tempj) * CLONG_rtz(d_Nx) + CLONG_rtz(tempi) + CLONG_rtz(tempk) * CLONG_rtz(d_Nyx);
	return local_ind;
}

DEVICE float voxelValue(const float t0, const float tc, const float L) {
	return (t0 - tc) * L;
}

// #ifdef SIDDON
// compute the distance that the ray traverses in the current voxel
DEVICE float compute_element(PTR_THR float* t0, PTR_THR float* tc, const float L, const float tu, const int u, PTR_THR int* temp_ijk, PTR_THR bool* pass) {
	*pass = (*t0 >= 0.f && *t0 <= 1.f) || (*tc >= 0.f && *tc <= 1.f);
	float local_ele = 0.f;
	if (*pass)
		local_ele = voxelValue(FMIN(*t0, 1.f), FMAX(*tc, 0.f), L);
	*temp_ijk += u;
	*tc = *t0;
	*t0 += tu;
	return local_ele;
}

// compute the initial voxel index (beginning of the ray)
DEVICE int voxel_index(const float pt, const float diff, const float d, const float apu) {
	return CINT_rtz((pt * diff - apu) / d);
}

DEVICE bool siddon_pre_loop_2D(const float b1, const float b2, const float diff1, const float diff2, const float max1, const float max2,
	const float d1, const float d2, const uint N1, const uint N2, PTR_THR int *temp1, PTR_THR int *temp2, PTR_THR float *t1u, PTR_THR float *t2u, PTR_THR uint *Np,
	const int TYYPPI, const float ys, const float xs, const float yd, const float xd, PTR_THR float *tc, PTR_THR int *u1, PTR_THR int *u2, PTR_THR float *t10, PTR_THR float *t20, PTR_THR bool *xy) {
	// If neither x- nor y-directions are perpendicular
// Correspond to the equations (9) and (10) from reference [2]
	const float apu_tx = b1 - xs;
	const float apu_ty = b2 - ys;
	*t10 = (apu_tx) / (diff1);
	*t20 = (apu_ty) / (diff2);
	const float txback = (max1 - xs) / (diff1);
	const float tyback = (max2 - ys) / (diff2);

	// Equations (5-8)
	const float txmin = FMIN(*t10, txback);
	const float txmax = FMAX(*t10, txback);
	const float tymin = FMIN(*t20, tyback);
	const float tymax = FMAX(*t20, tyback);

	// (3-4)
	*tc = FMAX(txmin, tymin);
	const float tmax = FMIN(txmax, tymax);
#ifdef ORTH
	if (*tc == *t10 || *tc == txback)
		*xy = true;
	else
		*xy = false;
#endif

	uint imin, imax, jmin, jmax;

		// If true, then the ray/LOR does not intersect the pixel space --> continue to the next LOR
		if (*tc >= tmax) {
			return true;
		}

		// (11-14)
		if (xs < xd)
			d_g_s_precomp(*tc, txmin, tmax, txmax, &imin, &imax, t10, u1, diff1, b1, d1, xs, N1);
		// (15-18)
		else
			s_g_d_precomp(*tc, txmin, tmax, txmax, &imin, &imax, t10, u1, diff1, b1, d1, xs, N1);

		//Same as above
		if (ys < yd)
			d_g_s_precomp(*tc, tymin, tmax, tymax, &jmin, &jmax, t20, u2, diff2, b2, d2, ys, N2);
		else
			s_g_d_precomp(*tc, tymin, tmax, tymax, &jmin, &jmax, t20, u2, diff2, b2, d2, ys, N2);

		*Np = imax + 1u + jmax + 1u - imin - jmin;

	// (2) and (19)
	const float pt = ((FMIN(*t10, *t20) + *tc) / 2.f);

	// (26)
	*temp1 = voxel_index(pt, diff1, d1, apu_tx);
	// (27)
	*temp2 = voxel_index(pt, diff2, d2, apu_ty);

	// (28)
	*t1u = d1 / FABS(diff1);
	*t2u = d2 / FABS(diff2);

	if (TYYPPI == 0) {
		if (*temp1 < 0 || *temp2 < 0 || *temp1 >= N1 || *temp2 >= N2)
			return true;
	}

	return false;
}

DEVICE bool siddon_pre_loop_3D(const FLOAT3 b, const FLOAT3 diff, const FLOAT3 max, const FLOAT3 dd, const uint3 N, PTR_THR int *tempi, PTR_THR int *tempj, PTR_THR int *tempk, 
    PTR_THR float *txu, PTR_THR float *tyu, PTR_THR float *tzu, PTR_THR uint *Np, const int TYYPPI, const FLOAT3 s, const FLOAT3 d, PTR_THR float *tc, PTR_THR int *i, PTR_THR int *j, PTR_THR int *k, PTR_THR float *tx0, 
	PTR_THR float *ty0, PTR_THR float *tz0, PTR_THR bool *xy, const int3 ii) {

	const float3 apuT = (float3)b - (float3)s;
	const float3 t0 = apuT / (float3)diff;
	const float3 tBack = DIVIDE3((float3)max - (float3)s, (float3)diff);

	const float3 tMin = FMIN(t0, tBack);
	const float3 tMax = FMAX(t0, tBack);

	*tc = FMAX(FMAX(tMin.x, tMin.z), tMin.y);
	const float tmax = FMIN(FMIN(tMax.x, tMax.z), tMax.y);
	*tx0 = t0.x;
	*ty0 = t0.y;
	*tz0 = t0.z;
#ifdef ORTH
	const FLOAT pituus = d.x - s.x;
	const FLOAT pituusY = d.y - s.y;
	const float angle = FABS(ACOS((pituus) / SQRT(pituus * pituus + pituusY * pituusY)));
	if ((angle < 0.785398f && angle > FLOAT_ZERO) || (angle > 2.35619f && angle < 3.92699f) || (angle > 5.497787f))
		*xy = true;
	else
		*xy = false;
#endif

	uint imin, imax, jmin, jmax, kmin, kmax;

		if (*tc >= tmax) {
			return true;
		}
		//float ax = 
		if (s.x < d.x)
			d_g_s_precomp(*tc, tMin.x, tmax, tMax.x, &imin, &imax, tx0, i, diff.x, b.x, dd.x, s.x, N.x);
		else
			s_g_d_precomp(*tc, tMin.x, tmax, tMax.x, &imin, &imax, tx0, i, diff.x, b.x, dd.x, s.x, N.x);

		if (s.y < d.y)
			d_g_s_precomp(*tc, tMin.y, tmax, tMax.y, &jmin, &jmax, ty0, j, diff.y, b.y, dd.y, s.y, N.y);
		else
			s_g_d_precomp(*tc, tMin.y, tmax, tMax.y, &jmin, &jmax, ty0, j, diff.y, b.y, dd.y, s.y, N.y);

		if (s.z < d.z)
			d_g_s_precomp(*tc, tMin.z, tmax, tMax.z, &kmin, &kmax, tz0, k, diff.z, b.z, dd.z, s.z, N.z);
		else
			s_g_d_precomp(*tc, tMin.z, tmax, tMax.z, &kmin, &kmax, tz0, k, diff.z, b.z, dd.z, s.z, N.z);

		*Np = (kmax - kmin + 1) + (jmax - jmin + 1) + (imax - imin + 1);

	const float pt = ((FMIN(FMIN(*tz0, *ty0), *tx0) + *tc) / 2.f);

	const float3 tempijkF = CLAMP3(FMAD3(pt, diff, -apuT) / dd, 0.f, CFLOAT3(N - 1));
	const int3 tempijk = CINT3_rtz(tempijkF);
	*tempi = tempijk.x;
	*tempj = tempijk.y;
	*tempk = tempijk.z;

	*txu = dd.x / (float)FABS(diff.x);
	*tyu = dd.y / (float)FABS(diff.y);
	*tzu = dd.z / (float)FABS(diff.z);

	return false;
}
#endif

#if defined(FP) && !defined(PROJ5)
DEVICE void forwardProjectAF(CLGLOBAL float* output, PTR_THR float *ax, size_t idx, const float temp, const int kk) {
    output[idx] += ax[kk]
#ifndef CT
	* temp
#endif
    ;
}
#endif

#ifdef HELICAL
DEVICE int rayArcIntersection(
    float3 s, 
    float3 v,
	float3 d, 
    float xc, float yc, float r,
    float zMin, float zMax,
    float thetaSpan, float thetaCenter, 
	float* theta, 
    float3* intersection
) {
	float fx = s.x - xc;
	float fy = s.y - yc;
    float a = v.x * v.x + v.y * v.y;
    float b = 2.f * (v.x * fx + v.y * fy);
    float c = fx * fx + fy * fy - r * r;
    float disc = b * b - 4.f * a * c;
    if (disc < 0) 
		return 0;

    for (int signv = -1; signv <= 1; signv += 2) {
        float t = (-b + (float)signv * SQRT(disc)) / (2.f*a);
        if (t < 0.f) 
			continue;
		float3 xyz = s + v * t;

		const float y = xyz.y - yc;
		const float x = xyz.x - xc;
		const float distX = xyz.x - d.x;
		const float distY = xyz.y - d.y;
		const float rr = 2.f * r * r;
		*theta = ACOS((rr - (distX * distX + distY * distY)) / rr);
        *theta *= -sign(ATAN2(y, x) - thetaCenter);
		if (*theta > thetaSpan / 2.f || *theta < 0)
			continue;
		float dist1 = xyz.x - s.x;
		float dist2 = xyz.y - s.y;
		float dist = SQRT(dist1 * dist1 + dist2 * dist2);
		float distV = SQRT(v.x * v.x + v.y * v.y);
		float angle = v.z/distV;
		xyz.z = angle * dist + d.z;
        if (xyz.z < zMin || xyz.z > zMax) 
			continue;

        *intersection = xyz;
		*theta += thetaSpan / 2.f;
        return 1;
    }
    return 0;
}

DEVICE void normalizeCurvedCoordinates(float3 xyz, 
	float zMin, float zMax,
    float thetaSpan,
	float theta, 
	float* u, float* v) {

    *u = theta / thetaSpan;

    *v = (xyz.z - zMin) / (zMax - zMin);
}
#endif

///////////////////////////////////////////////////////////////////////
///////////////// COMMON PRIOR AND ALGORITHM FUNCTIONS ////////////////
///////////////////////////////////////////////////////////////////////
// These were added to share as much code with fastPDHG as possible
#if defined(NLM_) || defined(FASTNLM) // START SHARED NLM

#ifndef NLTYPE
#define NLTYPE 0
#endif

// Read a single voxel of the input from an image/texture
#if defined(OPENCL)
CONSTANT sampler_t samplerNLMShared = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_CLAMP_TO_EDGE;
#endif
DEVICE float NLMTexRead(IMAGE3D uIn, const int x, const int y, const int z) {
#if defined(CUDA) || defined(HIP)
    return tex3D<float>(uIn, x, y, z);
#else
    return read_imagef(uIn, samplerNLMShared, (int4)(x, y, z, 0)).w;
#endif
}

// The NLM neighborhood is read either from a flat local-memory cache or straight from an image/texture
// The latter only when FASTNLMLOCAL is set as false
#if defined(NLM_) || defined(FASTNLMLOCAL)
#define NLM_USE_LOCAL 1
#endif

#ifdef NLM_USE_LOCAL
// z-extent of one work-group
// fastPDHG needs to take into account NVOXELS and similar also
#if defined(FASTNLM)
#define NLM_TILEZ_BASE FASTNLMTILEZ
#else
#define NLM_TILEZ_BASE LOCAL_SIZE3
#endif
#define NLM_TILEX (LOCAL_SIZE + SWINDOWX * 2 + PWINDOWX * 2)
#define NLM_TILEY (LOCAL_SIZE2 + SWINDOWY * 2 + PWINDOWY * 2)
#define NLM_TILEZ (NLM_TILEZ_BASE + SWINDOWZ * 2 + PWINDOWZ * 2)
#if defined(OPENCL)
#define NLM_LOCALQ __local
#else
#define NLM_LOCALQ
#endif
#define NLM_IN_T NLM_LOCALQ const float* CLRESTRICT
// Linear indexing with contiguous fetching
#define NLMFETCH(BUF, X, Y, Z) BUF[(X) + (Y) * NLM_TILEX + (Z) * NLM_TILEX * NLM_TILEY]
#else
#define NLM_IN_T IMAGE3D
#define NLMFETCH(BUF, X, Y, Z) NLMTexRead(BUF, X, Y, Z)
#endif

// Limit loop unrolling to 128 to avoid too large (machine) code and compilation times
#define NLM_ITER ((SWINDOWX * 2 + 1) * (SWINDOWY * 2 + 1) * (SWINDOWZ * 2 + 1) * (PWINDOWX * 2 + 1) * (PWINDOWY * 2 + 1) * (PWINDOWZ * 2 + 1))
#define NLM_UNROLL (NLM_ITER <= 128)

// Gradient of the NLM prior (and its variants) for a single voxel
DEVICE float NLMGradient(NLM_IN_T uIn,
    CONSTANT float* gaussian,
    const int x, const int y, const int z, const float h, const float epps
#if NLTYPE >= 3
    , const float gamma
#endif
#if NLTYPE == 6 // NLGGMRF
    , const float p, const float q, const float c
#endif
#ifdef NLMADAPTIVE
    , const float s
#endif
#ifdef NLMREF
    , NLM_IN_T uRefIn
#endif
) {
    float weight_sum = epps;
    float output = FLOAT_ZERO;
#if NLTYPE == 1
    float outputAla = epps;
#endif
    const float uj = NLMFETCH(uIn, x, y, z);
#if NLTYPE == 6
    // Precompute for NLGGMRF
    const float cpq = POWR(c, p - q);
#endif
#ifdef NLMADAPTIVE
    const float pSize = CFLOAT((PWINDOWX * 2 + 1) * (PWINDOWY * 2 + 1) * (PWINDOWZ * 2 + 1));
    float hh = FLOAT_ZERO;
#endif
    // The patch around the voxel itself is the same for every offset of the search window, so keep it in registers
    // (private memory cache) rather than re-reading it
    float pCache[(PWINDOWX * 2 + 1) * (PWINDOWY * 2 + 1) * (PWINDOWZ * 2 + 1)];
    {
        int pInd = 0;
#pragma unroll
        for (int pz = -PWINDOWZ; pz <= PWINDOWZ; pz++) {
#pragma unroll
            for (int py = -PWINDOWY; py <= PWINDOWY; py++) {
#pragma unroll
                for (int px = -PWINDOWX; px <= PWINDOWX; px++) {
#ifdef NLMREF
                    pCache[pInd++] = NLMFETCH(uRefIn, x + px, y + py, z + pz);
#else
                    pCache[pInd++] = NLMFETCH(uIn, x + px, y + py, z + pz);
#endif
                }
            }
        }
    }
#if NLM_UNROLL
#pragma unroll
#endif
    for (int i = -SWINDOWX; i <= SWINDOWX; i++) {
#if NLM_UNROLL
#pragma unroll
#endif
        for (int j = -SWINDOWY; j <= SWINDOWY; j++) {
#if NLM_UNROLL
#pragma unroll
#endif
            for (int k = -SWINDOWZ; k <= SWINDOWZ; k++) {
                if (i == 0 && j == 0 && k == 0)
                    continue;
                float weight = FLOAT_ZERO;
                float distance = FLOAT_ZERO;
                int pInd = 0;
#pragma unroll
                for (int pz = -PWINDOWZ; pz <= PWINDOWZ; pz++) {
#pragma unroll
                    for (int py = -PWINDOWY; py <= PWINDOWY; py++) {
                        int dim_g = (pz + PWINDOWZ) * (PWINDOWX * 2 + 1) * (PWINDOWY * 2 + 1) + (py + PWINDOWY) * (PWINDOWX * 2 + 1);
#pragma unroll
                        for (int px = -PWINDOWX; px <= PWINDOWX; px++) {
                            const float gg = gaussian[dim_g++];
#ifdef NLMREF
                            const float Pk = NLMFETCH(uRefIn, x + i + px, y + j + py, z + k + pz);
#else
                            const float Pk = NLMFETCH(uIn, x + i + px, y + j + py, z + k + pz);
#endif
                            const float PP = pCache[pInd++] - Pk;
                            distance += gg * PP * PP;
                        }
                    }
                }
#ifdef NLMADAPTIVE
                hh = distance / pSize;
                weight = EXP(-distance / (hh * h + s));
#else
                weight = EXP(-distance / h);
#endif
                weight_sum += weight;
                const float uk = NLMFETCH(uIn, x + i, y + j, z + k);
                // Different NLM regularization methods
                // NLTYPE 0 = MRF NLM
                // NLTYPE 1 = NLTV
                // NLTYPE 2 = NLM filtered (i.e. similar to MRP)
                // NLTYPE 3 = NLRD
                // NLTYPE 4 = NL Lange
                // NLTYPE 5 = NLM filtered with Lange
                // NLTYPE 6 = NLGGMRF
                // NLTYPE 7 = ?
#if NLTYPE == 2 || NLTYPE == 5 // START NLM NLTYPE
                // NLMRP
                output += weight * uk;
#elif NLTYPE == 0
                output += (weight * (uj - uk));
#elif NLTYPE == 3
                // NLRD
                const float u = (uj - uk);
#ifndef USEMAD // START FMAD
                const float divPow = (uj + uk + gamma * fabs(u) + epps);
                output += weight * u * (gamma * fabs(u) + uj + 3.f * uk + epps * epps) / (divPow * divPow);
#else
                const float divPow = FMAD(gamma, fabs(u), uj + uk + epps);
                output += weight * u * (FMAD(gamma, fabs(u), uj + 3.f * uk + epps * epps)) / (divPow * divPow);
#endif // END FMAD
#elif NLTYPE == 4
                // Lange
                const float u = (uj - uk);
                const float uabs = sign(u);
                output += weight * (uabs - uabs / (fabs(u) / gamma + FLOAT_ONE));
#elif NLTYPE == 6
                // NLGGMRF
                const float delta = uj - uk;
                const float dcpq = POWR(fabs(delta / c), p - q);
                const float deltapqc = FLOAT_ONE + dcpq;
                output += weight * (POWR(fabs(delta), p - FLOAT_ONE) / deltapqc) * (p - gamma * ((dcpq * cpq) / deltapqc)) * sign(delta);
#elif NLTYPE == 7
                const float u = (uk - uj);
                const float apu = (u * u + gamma * gamma);
                output += ((FLOAT_TWO * u * u * u) / (apu * apu) - FLOAT_TWO * (u / apu));
#else
                //NLTV
                const float apuU = uj - uk;
                output += (weight * apuU);
                outputAla += weight * apuU * apuU;
#endif // END NLM NLTYPE
            }
        }
    }
    weight_sum = FLOAT_ONE / weight_sum;
    output *= weight_sum;
#if NLTYPE == 2 // START NLM NLTYPE
    output = uj - output;
#elif NLTYPE == 5
    // Lange with NLMRP
    output = uj - output;
    const float uabs = sign(output);
    output = (uabs - uabs / (fabs(output) / gamma + FLOAT_ONE));
#elif NLTYPE == 1
#ifndef USEMAD // START FMAD
    output /= SQRT(outputAla * weight_sum + epps);
#else
    output /= SQRT(FMAD(outputAla, weight_sum, epps));
#endif // END FMAD
#endif // END NLM NLTYPE
    return output;
}

#endif // END SHARED NLM

// Shared search-window prior gradient (RDP, GGMRF, hyperbolic)
#if defined(RDP) || defined(GGMRF) || defined(HYPER) || defined(FASTRDP) || defined(FASTGGMRF) || defined(FASTHYPER) || defined(TVGRAD) || defined(FASTTV) // START SHARED SW PRIOR

// Reference image support
#if defined(RDPREF) && !defined(PRIORREF)
#define PRIORREF
#endif

// Sampler shared by the search-window priors (RDP/GGMRF/hyperbolic)
#if defined(OPENCL)
CONSTANT sampler_t samplerSWPrior = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_CLAMP_TO_EDGE;
#endif
DEVICE float SWTexRead(IMAGE3D uIn, const int x, const int y, const int z) {
#if defined(CUDA) || defined(HIP)
    return tex3D<float>(uIn, x, y, z);
#else
    return read_imagef(uIn, samplerSWPrior, (int4)(x, y, z, 0)).w;
#endif
}
DEVICE float SWBufRead(const CLGLOBAL float* CLRESTRICT buf, const int x, const int y, const int z, const int3 N) {
    if (x < 0 || y < 0 || z < 0 || x >= N.x || y >= N.y || z >= N.z)
        return FLOAT_ZERO;
    return buf[(x) + (y) * N.x + (z) * N.x * N.y];
}


#if defined(RDP) || defined(FASTRDP) // START SHARED RDP GRADIENT
// Both image/texture and pure buffer support
#if defined(FASTRDP) || (defined(RDP) && defined(USEIMAGES))
#define RDPFETCH(U, X, Y, Z) SWTexRead(U, X, Y, Z)
#define RDP_IN_T IMAGE3D
#else
#define RDPFETCH(U, X, Y, Z) SWBufRead(U, X, Y, Z, rdpN)
#define RDP_IN_T const CLGLOBAL float* CLRESTRICT
#define RDPNEEDN
#endif

// Standard RDP which only takes the (6) neighboring voxels into account
DEVICE float RDPGradientNorm(RDP_IN_T u, const int x, const int y, const int z, const float gamma, const float epps
#ifdef RDPNEEDN
    , const int3 rdpN
#endif
) {
    const float uj = RDPFETCH(u, x, y, z);
    const float2 ux = MFLOAT2(RDPFETCH(u, x + 1, y, z), RDPFETCH(u, x - 1, y, z));
    const float2 uy = MFLOAT2(RDPFETCH(u, x, y + 1, z), RDPFETCH(u, x, y - 1, z));
    const float2 uz = MFLOAT2(RDPFETCH(u, x, y, z + 1), RDPFETCH(u, x, y, z - 1));
    const float2 uj_ux = uj - ux;
    const float2 uj_uy = uj - uy;
    const float2 uj_uz = uj - uz;
#ifndef USEMAD
    const float2 divPow2X = (uj + ux + gamma * fabs(uj_ux));
    const float2 divPow2Y = (uj + uy + gamma * fabs(uj_uy));
    const float2 divPow2Z = (uj + uz + gamma * fabs(uj_uz));
    float2 output = uj_ux * (gamma * fabs(uj_ux) + uj + 3.f * ux + epps * epps) / (divPow2X * divPow2X + epps)
        + uj_uy * (gamma * fabs(uj_uy) + uj + 3.f * uy + epps * epps) / (divPow2Y * divPow2Y + epps)
        + uj_uz * (gamma * fabs(uj_uz) + uj + 3.f * uz + epps * epps) / (divPow2Z * divPow2Z + epps);
#else
    const float2 divPow2X = FMAD2(gamma, fabs(uj_ux), uj + ux);
    const float2 divPow2Y = FMAD2(gamma, fabs(uj_uy), uj + uy);
    const float2 divPow2Z = FMAD2(gamma, fabs(uj_uz), uj + uz);
    float2 output = uj_ux * FMAD2(gamma, fabs(uj_ux), uj + 3.f * ux + epps * epps) / (divPow2X * divPow2X + epps)
        + uj_uy * FMAD2(gamma, fabs(uj_uy), uj + 3.f * uy + epps * epps) / (divPow2Y * divPow2Y + epps)
        + uj_uz * FMAD2(gamma, fabs(uj_uz), uj + 3.f * uz + epps * epps) / (divPow2Z * divPow2Z + epps);
#endif
    if (isnan(output.x))
		output.x = FLOAT_ZERO;
    if (isnan(output.y))
		output.y = FLOAT_ZERO;
    return output.x + output.y;
}
#endif // END SHARED RDP GRADIENT

// RDP "with corners" uses the same predefined neighborhood as GGMRF and hyperbolic
#if defined(RDPCORNERS) || defined(GGMRF) || defined(HYPER) || defined(FASTRDPCORNERS) || defined(FASTGGMRF) || defined(FASTHYPER)

// Whether local memory caching is used
#if defined(RDPCORNERS) || defined(GGMRF) || defined(HYPER) || defined(FASTSWLOCAL)
#define SW_USE_LOCAL 1
#endif
#ifdef SW_USE_LOCAL
// The axial (z) cache is different depending on whether fastPDHG is used or not
// fastPDHG combines this with backprojection which can compute multiple slices in the same thread
#if defined(FASTSWLOCAL)
#define SW_TILEZ_BASE FASTNLMTILEZ
#else
#define SW_TILEZ_BASE LOCAL_SIZE3
#endif
#define SW_TILEX (LOCAL_SIZE + SWINDOWX * 2)
#define SW_TILEY (LOCAL_SIZE2 + SWINDOWY * 2)
#define SW_TILEZ (SW_TILEZ_BASE + SWINDOWZ * 2)
#if defined(OPENCL)
#define SW_LOCALQ __local
#else
#define SW_LOCALQ
#endif
#define SW_IN_T SW_LOCALQ const float* CLRESTRICT
#define SWFETCH(BUF, X, Y, Z) BUF[(X) + (Y) * SW_TILEX + (Z) * SW_TILEX * SW_TILEY]
#else
#define SW_IN_T IMAGE3D
#define SWFETCH(BUF, X, Y, Z) SWTexRead(BUF, X, Y, Z)
#endif

// RDP, hyperbolic and GGMRF gradients
// All compute it in the specified neighborhood (Ndx/y/z)
// Common function for both fastPDHG and "regular" cases
// Hyperbolic uses different ordering with the weights (z --> y --> x) than GGMRF and RDP
DEVICE float priorGradientSW(SW_IN_T uIn, CONSTANT float* weight,
    const int x, const int y, const int z, const float epps
#if SWPRIORTYPE == 0
    , const float gamma
#ifdef PRIORREF
    , SW_IN_T uRefIn
#endif
#endif
#if SWPRIORTYPE == 1   // GGMRF
    , const float p, const float q, const float c, const float pqc
#endif
#if SWPRIORTYPE == 2   // hyperbolic
    , const float sigma
#endif
) {
    float output = FLOAT_ZERO;
    const float uj = SWFETCH(uIn, x, y, z);
#if SWPRIORTYPE == 0 && defined(PRIORREF)
    const float kj = SWFETCH(uRefIn, x, y, z);
#endif
#if SWPRIORTYPE == 1
    const float cpq = POWR(c, p - q);
#endif
#if SWPRIORTYPE == 2
    const float invSigma = FLOAT_ONE / sigma;
#endif
    int uu = 0;
#if SWPRIORTYPE == 2 // START HYPER
    for (int k = -SWINDOWZ; k <= SWINDOWZ; k++) {
        for (int j = -SWINDOWY; j <= SWINDOWY; j++) {
            for (int i = -SWINDOWX; i <= SWINDOWX; i++) {
                if (i == 0 && j == 0 && k == 0)
                    continue;
                const float uk = SWFETCH(uIn, x + i, y + j, z + k);
                const float ux = (uj - uk) * invSigma;
                output += ux * invSigma * RSQRT(FMAD(ux, ux, FLOAT_ONE)) * weight[uu];
                uu++;
            }
        }
    }
#else // START RDP/GGMRF
    for (int i = -SWINDOWX; i <= SWINDOWX; i++) {
        for (int j = -SWINDOWY; j <= SWINDOWY; j++) {
            for (int k = -SWINDOWZ; k <= SWINDOWZ; k++) {
                if (i == 0 && j == 0 && k == 0)
                    continue;
                const float uk = SWFETCH(uIn, x + i, y + j, z + k);
#if SWPRIORTYPE == 0 // RDP
                const float delta = uj - uk;
                const float divPow2 = FMAD(gamma, fabs(delta), uj + uk);
#ifdef PRIORREF
                const float kk = SWFETCH(uRefIn, x + i, y + j, z + k);
                output += weight[uu] * SQRT(kk * kj) * delta * (gamma * fabs(delta) + uj + 3.f * uk + epps * epps) / (divPow2 * divPow2 + epps);
#else
                output += weight[uu] * delta * (gamma * fabs(delta) + uj + 3.f * uk + epps * epps) / (divPow2 * divPow2 + epps);
#endif
#endif
#if SWPRIORTYPE == 1 // GGMRF
                const float delta = uj - uk;
                if (delta != FLOAT_ZERO) {
                    const float dcpq = POWR(fabs(delta / c), p - q);
                    const float deltapqc = FLOAT_ONE + dcpq;
                    output += weight[uu] * (POWR(fabs(delta), p - FLOAT_ONE) / deltapqc) * (p - pqc * ((dcpq * cpq) / deltapqc)) * sign(delta);
                }
#endif
                uu++;
            }
        }
    }
#endif // END HYPER / RDP+GGMRF loop-order dispatch
#if SWPRIORTYPE == 0
    if (isnan(output))
        output = FLOAT_ZERO;
#endif
    return output;
}

#endif // END SW window (RDPCORNERS/GGMRF/HYPER/FASTRDPCORNERS/FASTGGMRF/FASTHYPER)

#endif // END SHARED SW PRIOR

// Common TV functions
// Anatomical reference image cases are not supported for fastPDHG
#if defined(TVGRAD) || defined(FASTTV)

// Image/texture version and buffer version
#if defined(FASTTV) || (defined(TVGRAD) && defined(USEIMAGES))
#define TVFETCH(U, X, Y, Z) SWTexRead(U, X, Y, Z)
#define TV_IN_T IMAGE3D
#else
#define TVFETCH(U, X, Y, Z) SWBufRead(U, X, Y, Z, tvN)
#define TV_IN_T const CLGLOBAL float* CLRESTRICT
#endif
// Load the volume dimensions if needed (buffer version and reference image cases)
#if !defined(USEIMAGES) || defined(ANATOMICAL1) || defined(ANATOMICAL2) || defined(ANATOMICAL3)
#define TVNEEDN
#endif

// Moved from auxKernels.cl
DEVICE float sqrtVal(const float3 input, const float epps
#ifdef TVW1
    , const float3 w
#endif
) {
#ifdef TVW1
#ifdef USEMAD
    return SQRT(FMAD(w.x, input.x * input.x, FMAD(w.y, input.y * input.y, FMAD(w.z, input.z * input.z, epps))));
#else
    return SQRT(w.x * input.x * input.x + w.y * input.y * input.y + w.z * input.z * input.z + epps);
#endif
#else
#ifdef USEMAD
    return SQRT(FMAD(input.x, input.x, FMAD(input.y, input.y, FMAD(input.z, input.z, epps))));
#else
    return SQRT(input.x * input.x + input.y * input.y + input.z * input.z + epps);
#endif
#endif
}

// Gradient of the TV prior for one voxel, covering every variant
DEVICE float TVGradient(TV_IN_T u, const int x, const int y, const int z, const float sigma, const float epps
#ifdef TVNEEDN
    , const int3 tvN
#endif
#if defined(ANATOMICAL2) || defined(ANATOMICAL3)
    , const float C
#endif
#if defined(ANATOMICAL1) || defined(ANATOMICAL2) || defined(ANATOMICAL3)
    , CLGLOBAL float* CLRESTRICT S
#endif
) {
    const float uijk = TVFETCH(u, x, y, z);
#if defined(SATV) // START JPTV || SATV
    float2 ux = MFLOAT2(TVFETCH(u, x + 1, y, z), TVFETCH(u, x - 1, y, z));
    float2 uy = MFLOAT2(TVFETCH(u, x, y + 1, z), TVFETCH(u, x, y - 1, z));
    float2 uz = MFLOAT2(TVFETCH(u, x, y, z + 1), TVFETCH(u, x, y, z - 1));
    ux = uijk - ux;
    uy = uijk - uy;
    uz = uijk - uz;
    const float2 uabsx = ux / (fabs(ux) + epps);
    const float2 uabsy = uy / (fabs(uy) + epps);
    const float2 uabsz = uz / (fabs(uz) + epps);
    float2 output = uabsx - uabsx / (fabs(ux) / sigma + FLOAT_ONE) + uabsy - uabsy / (fabs(uy) / sigma + FLOAT_ONE) + uabsz - uabsz / (fabs(uz) / sigma + FLOAT_ONE);
    return output.x + output.y;
#else
    const float3 uijkP = MFLOAT3(TVFETCH(u, x + 1, y, z), TVFETCH(u, x, y + 1, z), TVFETCH(u, x, y, z + 1));
    const float3 uijkM = MFLOAT3(TVFETCH(u, x - 1, y, z), TVFETCH(u, x, y - 1, z), TVFETCH(u, x, y, z - 1));
    const float2 ui = MFLOAT2(TVFETCH(u, x - 1, y + 1, z), TVFETCH(u, x - 1, y, z + 1));
    const float2 uj = MFLOAT2(TVFETCH(u, x + 1, y - 1, z), TVFETCH(u, x, y - 1, z + 1));
    const float2 uk = MFLOAT2(TVFETCH(u, x + 1, y, z - 1), TVFETCH(u, x, y + 1, z - 1));
    const float3 u1 = MFLOAT3(uijk - uijkM.x, ui.x - uijkM.x, ui.y - uijkM.x);
    const float3 u2 = MFLOAT3(uj.x - uijkM.y, uijk - uijkM.y, uj.y - uijkM.y);
    const float3 u3 = MFLOAT3(uk.x - uijkM.z, uk.y - uijkM.z, uijk - uijkM.z);
#ifdef TVW1 // START TVW1
    const float3 u4 = uijkP - uijk;
    float3 w4 = (u4) / sigma;
    w4 = EXP3(-w4 * w4);
    const float pvalijk = sqrtVal(u4, epps, w4);
    float3 w1 = (u1) / sigma;
    w1 = EXP3(-w1 * w1);
    float3 w2 = (u2) / sigma;
    w2 = EXP3(-w2 * w2);
    float3 w3 = (u3) / sigma;
    w3 = EXP3(-w3 * w3);
#ifdef USEMAD
    return -(FMAD(w4.x, u4.x, FMAD(w4.y, u4.y, w4.z * u4.z))) / pvalijk + (w1.x * (uijk - uijkM.x)) / sqrtVal(u1, epps, w1) + (w2.y * (uijk - uijkM.y)) / sqrtVal(u2, epps, w2) + (w3.y * (uijk - uijkM.z)) / sqrtVal(u3, epps, w3);
#else
    return -(w4.x * u4.x + w4.y * u4.y + w4.z * u4.z) / pvalijk + (w1.x * (uijk - uijkM.x)) / sqrtVal(u1, epps, w1) + (w2.y * (uijk - uijkM.y)) / sqrtVal(u2, epps, w2) + (w3.y * (uijk - uijkM.z)) / sqrtVal(u3, epps, w3);
#endif
#else
#ifdef ANATOMICAL1 // TV type 1
    const int NN = tvN.x * tvN.y * tvN.z;
    const int n = x + y * tvN.x + z * tvN.x * tvN.y;
#if defined(CUDA) || defined(HIP)
    float s[9];
#else
    __private float s[9];
#endif
    for (int kk = 0; kk < 9; kk++)
        s[kk] = S[n + NN * kk];
    const float3 val = uijkP - uijk;
    const float pvalijk = SQRT(val.x * val.x * s[0] + val.y * val.y * s[4] + val.z * val.z * s[8] + s[1] * (val.x) * (val.y) + s[3] * (val.x) * (val.y) + s[2] * (val.x) * (val.z) + s[6] * (val.x) * (val.z) +
        s[5] * (val.y) * (val.z) + s[7] * (val.y) * (val.z) + epps);
    const float pvalijkX = SQRT(u1.x * u1.x * s[0] + u1.y * u1.y * s[4] + u1.z * u1.z * s[8] + s[1] * (u1.x) * (u1.y) + s[3] * (u1.x) * (u1.y) + s[2] * (u1.x) * (u1.z) + s[6] * (u1.x) * (u1.z) +
        s[5] * (u1.y) * (u1.z) + s[7] * (u1.y) * (u1.z) + epps);
    const float pvalijkY = SQRT(u2.x * u2.x * s[0] + u2.y * u2.y * s[4] + u2.z * u2.z * s[8] + s[1] * (u2.x) * (u2.y) + s[3] * (u2.x) * (u2.y) + s[2] * (u2.x) * (u2.z) + s[6] * (u2.x) * (u2.z) +
        s[5] * (u2.y) * (u2.z) + s[7] * (u2.y) * (u2.z) + epps);
    const float pvalijkZ = SQRT(u3.x * u3.x * s[0] + u3.y * u3.y * s[4] + u3.z * u3.z * s[8] + s[1] * (u3.x) * (u3.y) + s[3] * (u3.x) * (u3.y) + s[2] * (u3.x) * (u3.z) + s[6] * (u3.x) * (u3.z) +
        s[5] * (u3.y) * (u3.z) + s[7] * (u3.y) * (u3.z) + epps);
    const float dx = s[0] * (FLOAT_TWO * (uijk - uijkM.x)) + s[3] * u1.y + s[2] * u1.z + s[6] * u1.z + s[1] * u1.y;
    const float dy = s[4] * (FLOAT_TWO * (uijk - uijkM.y)) + s[5] * u2.z + s[3] * u2.x + s[1] * u2.x + s[7] * u2.z;
    const float dz = s[8] * (FLOAT_TWO * (uijk - uijkM.z)) + s[6] * u3.x + s[5] * u3.y + s[7] * u3.y + s[2] * u3.x;
    const float d = s[1] * val.x + s[2] * val.x + s[3] * val.x + s[6] * val.x + s[1] * val.y + s[3] * val.y + s[5] * val.y + s[7] * val.y + s[2] * val.z + s[5] * val.z + s[6] * val.z + s[7] * val.z + s[0] * FLOAT_TWO * val.x + s[4] * FLOAT_TWO * val.y + s[8] * FLOAT_TWO * val.z;
    return FLOAT_HALF * (d / pvalijk + dx / pvalijkX + dy / pvalijkY + dz / pvalijkZ);
#elif defined(ANATOMICAL2) // TV type 2
    const float3 uijkR = MFLOAT3(SWBufRead(S, x + 1, y, z, tvN), SWBufRead(S, x, y + 1, z, tvN), SWBufRead(S, x, y, z + 1, tvN));
    const float3 apuS = (uijkR - SWBufRead(S, x, y, z, tvN));
    const float3 apu = uijkP - uijk;
    const float pvalijk = SQRT(apu.x * apu.x + apu.y * apu.y + apu.z * apu.z + C * (apuS.x * apuS.x + apuS.y * apuS.y + apuS.z * apuS.z) + epps);
    return (3.f * uijk - uijkP.x - uijkP.y - uijkP.z) / pvalijk + (uijk - uijkM.x) / sqrtVal(u1, epps) + (uijk - uijkM.y) / sqrtVal(u2, epps) + (uijk - uijkM.z) / sqrtVal(u3, epps) + 1e-7f;
#elif defined(ANATOMICAL3) // APLS
    const float3 uijkR = MFLOAT3(SWBufRead(S, x + 1, y, z, tvN), SWBufRead(S, x, y + 1, z, tvN), SWBufRead(S, x, y, z + 1, tvN));
    float3 epsilon = (uijkR - SWBufRead(S, x, y, z, tvN));
    epsilon = epsilon / SQRT(epsilon.x * epsilon.x + epsilon.y * epsilon.y + epsilon.z * epsilon.z + C * C);
    const float3 apu = uijkP - uijk;
    const float apuR = uijkR.x * apu.x + uijkR.y * apu.y + uijkR.z * apu.z;
    const float pvalijk = SQRT(apu.x * apu.x + apu.y * apu.y + apu.z * apu.z - apuR * apuR + epps);
    float apuRXYZ = uijkR.x * u1.x + uijkR.y * u1.y + uijkR.z * u1.z;
    const float pvalijkX = SQRT(u1.x * u1.x + u1.y * u1.y + u1.z * u1.z + apuRXYZ * apuRXYZ + epps);
    apuRXYZ = uijkR.x * u2.x + uijkR.y * u2.y + uijkR.z * u2.z;
    const float pvalijkY = SQRT(u2.x * u2.x + u2.y * u2.y + u2.z * u2.z + apuRXYZ * apuRXYZ + epps);
    apuRXYZ = uijkR.x * u3.x + uijkR.y * u3.y + uijkR.z * u3.z;
    const float pvalijkZ = SQRT(u3.x * u3.x + u3.y * u3.y + u3.z * u3.z + apuRXYZ * apuRXYZ + epps);
    return FLOAT_HALF * ((6.f * uijk - FLOAT_TWO * uijkP.x - FLOAT_TWO * uijkP.y - FLOAT_TWO * uijkP.z + FLOAT_TWO * (epsilon.x*(uijk - uijkP.x) + epsilon.y*(uijk - uijkP.y) + epsilon.z*(uijk - uijkP.z)) * (epsilon.x + epsilon.y + epsilon.z)) / pvalijk +
        FLOAT_TWO * (u1.x - epsilon.x * (epsilon.x * u1.x + epsilon.y * u1.y + epsilon.z * u1.z)) / pvalijkX + FLOAT_TWO * (u2.y - epsilon.y * (epsilon.x * u2.x + epsilon.y * u2.y + epsilon.z * u2.z)) / pvalijkY +
        FLOAT_TWO * (u3.z - epsilon.z * (epsilon.x * u3.x + epsilon.y * u3.y + epsilon.z * u3.z))/ pvalijkZ + 1e-7f);
#else // Non-reference image TV (standard / JPTV --> identical)
    const float pvalijk = sqrtVal(uijkP - uijk, epps);
    return (3.f * uijk - uijkP.x - uijkP.y - uijkP.z) / pvalijk + (uijk - uijkM.x) / sqrtVal(u1, epps) + (uijk - uijkM.y) / sqrtVal(u2, epps) + (uijk - uijkM.z) / sqrtVal(u3, epps) + 1e-7f;
#endif
#endif // END TVW1
#endif // END JPTV || SATV
}

#endif // END TVGRAD || FASTTV

// PDHG subset image estimate update
#if defined(PDHG) || defined(FASTPDHG)
DEVICE float PDHGSubsetPrimal(const float imEst, const float rhs, const float tau, const float epps, const uchar enforcePositivity) {
    float imNew = imEst - tau * rhs;
    if (enforcePositivity != 0u)
        imNew = FMAX(epps, imNew);
    return imNew;
}
#endif

// PKMA/MBSREM/BSREM image estimate update
#if defined(PKMA) || defined(MBSREM) || defined(BSREM)
DEVICE float PoissonUpdateVoxel(const float imOld, const float rhs, const float lambda, const float epps,
	const float alpha, const uchar enforcePositivity) {
#ifdef PKMA
	float imApu = imOld - lambda * rhs;
#elif defined(MBSREM)
	float imApu = imOld + lambda * rhs;
#elif defined(BSREM)
	float imApu = imOld + lambda * rhs * imOld;
#endif
	if (enforcePositivity)
		imApu = fmax(epps, imApu);
#ifdef PKMA
	return (FLOAT_ONE - alpha) * imOld + alpha * imApu;
#elif defined(MBSREM)
	if (imApu >= alpha)
		imApu = alpha - epps;
	return imApu;
#elif defined(BSREM)
	return imApu;
#endif
}
#endif

// Set defaults if not defined
#ifdef FASTPDHG
#ifndef FASTPRECOND
// 0 = no image-based preconditioner, 1 = diagonal normalization (type 0), 2 = EM (type 1)
#define FASTPRECOND 0
#endif
#ifndef FASTALG
// 0 = PDHG, 1 = Poisson (PKMA / MBSREM / BSREM)
#define FASTALG 0
#endif

#if defined(FASTNLM) && defined(FASTNLMLOCAL)
// Fill the flat local-memory NLM neighborhood cache for the whole work-group
// Every work-item must call this (before the bounds check) so the following barrier is uniform
// iz is the work-item's global z
DEVICE void fastFillNLMCache(NLM_LOCALQ float* lCacheF,
#ifdef NLMREF
    NLM_LOCALQ float* lCacheRefF,
#endif
    IMAGE3D d_imNLM,
#ifdef NLMREF
    IMAGE3D d_urefNLM,
#endif
    const int iz, const int fastPriorZOffset
) {
    const int baseX = CINT(GRID0) * CINT(LSIZE0) - (SWINDOWX + PWINDOWX);
    const int baseY = CINT(GRID1) * CINT(LSIZE1) - (SWINDOWY + PWINDOWY);
    const int baseZ = iz + fastPriorZOffset - (SWINDOWZ + PWINDOWZ);
    const int nWorkItems = CINT(LSIZE0) * CINT(LSIZE1);
    const int lIndex = CINT(LID0) + CINT(LID1) * CINT(LSIZE0);
    for (int n = lIndex; n < NLM_TILEX * NLM_TILEY * NLM_TILEZ; n += nWorkItems) {
        const int cx = n % NLM_TILEX;
        const int cy = (n / NLM_TILEX) % NLM_TILEY;
        const int cz = n / (NLM_TILEX * NLM_TILEY);
        lCacheF[n] = NLMTexRead(d_imNLM, baseX + cx, baseY + cy, baseZ + cz);
#ifdef NLMREF
        lCacheRefF[n] = NLMTexRead(d_urefNLM, baseX + cx, baseY + cy, baseZ + cz);
#endif
    }
}
#endif // FASTNLM && FASTNLMLOCAL

#if defined(FASTSWLOCAL)
// Fill the flat local-memory neighborhood
// Affects RDP (with corners), GGMRF and hyperbolic
DEVICE void fastFillSWCache(SW_LOCALQ float* lCacheSW, IMAGE3D d_imNLM, const int iz, const int fastPriorZOffset) {
    const int baseX = CINT(GRID0) * CINT(LSIZE0) - SWINDOWX;
    const int baseY = CINT(GRID1) * CINT(LSIZE1) - SWINDOWY;
    const int baseZ = iz + fastPriorZOffset - SWINDOWZ;
    const int nWorkItems = CINT(LSIZE0) * CINT(LSIZE1);
    const int lIndex = CINT(LID0) + CINT(LID1) * CINT(LSIZE0);
    for (int n = lIndex; n < SW_TILEX * SW_TILEY * SW_TILEZ; n += nWorkItems) {
        const int cx = n % SW_TILEX;
        const int cy = (n / SW_TILEX) % SW_TILEY;
        const int cz = n / (SW_TILEX * SW_TILEY);
        lCacheSW[n] = SWTexRead(d_imNLM, baseX + cx, baseY + cy, baseZ + cz);
    }
}
#endif

// The image-domain part of fastPDHG
// Once backprojection has been computed this part is run
// Supports only PDHG (and its variants), PKMA/MBSREM/BSREM, RDP/GGMRF/NLM(and its variants)/hyperbolic/TV(gradient)
// As well as the first two image-based preconditioners
// Used automatically if applicable
DEVICE void fastPDHGUpdate(
    const float* temp, const float* wSum,
    const int3 i, size_t idx, const uint3 d_N, const int nVoxels, const uchar no_norm,
    // ii = volume number (0 = main volume)
	// The spatial prior is only applied to the main volume
    const int ii,
    // largeDim offsets
    const LONG fastImOffset, const int fastPriorZOffset,
#if FASTALG == 0
    CLGLOBAL float* CLRESTRICT d_U,
#endif
    CLGLOBAL float* CLRESTRICT d_OSEM, CLGLOBAL CAST* CLRESTRICT d_Summ,
#if FASTPRECOND > 0
    const CLGLOBAL float* CLRESTRICT d_precond,
#endif
#ifdef FASTNLM
    CONSTANT float* d_gaussian,
#ifdef FASTNLMLOCAL
    NLM_LOCALQ const float* CLRESTRICT lCacheF,
#ifdef NLMREF
    NLM_LOCALQ const float* CLRESTRICT lCacheRefF,
#endif
#else
    IMAGE3D d_imNLM,
#ifdef NLMREF
    IMAGE3D d_urefNLM,
#endif
#endif
    const float fastNLMh,
#if NLTYPE >= 3
    const float fastNLMgamma,
#endif
#if NLTYPE == 6
    const float fastNLMp, const float fastNLMq, const float fastNLMc,
#endif
#ifdef NLMADAPTIVE
    const float fastNLMs,
#endif
#endif // FASTNLM
#ifdef FASTSWLOCAL
    SW_LOCALQ const float* CLRESTRICT lCacheSW,
#endif
#ifdef FASTRDP
    IMAGE3D d_imNLM,
#ifdef FASTRDPCORNERS
    CONSTANT float* d_swWeight,
#endif
    const float fastRDPgamma,
#ifdef PRIORREF
    IMAGE3D d_urefNLM,
#endif
#endif // FASTRDP
#ifdef FASTGGMRF
    IMAGE3D d_imNLM,
    CONSTANT float* d_swWeight,
    const float fastGGMRFp, const float fastGGMRFq, const float fastGGMRFc, const float fastGGMRFpqc,
#endif // FASTGGMRF
#ifdef FASTHYPER
    IMAGE3D d_imNLM,
    CONSTANT float* d_swWeight,
    const float fastHYPERsigma,
#endif // FASTHYPER
#ifdef FASTTV
    // Anatomical reference image versions are not supported
    IMAGE3D d_imNLM,
    const float fastTVsigma,
#endif // FASTTV
#if FASTALG == 0
    const float fastTheta, const float fastTau,
#else
    const float fastLambda, const float fastAlpha,
#endif
    const float fastBeta,
    const float fastEpps, const uchar fastPositivity) {
    for (int zz = 0; zz < nVoxels; zz++) {
        const uint ind = i.z + zz;
        if (ind >= d_N.z)
            break;
        const float bp = temp[zz];
        const float curEst = d_OSEM[idx + fastImOffset];
#if FASTALG == 0
        // PDHG with subsets
        const float uNew = d_U[idx] + bp;
        d_U[idx] = uNew;
        float rhs = uNew + fastTheta * bp;
#else
		// Poisson
        float rhs = bp;
#endif
        // Compute the regularization if applicable
#if defined(FASTNLM) || defined(FASTRDP) || defined(FASTGGMRF) || defined(FASTHYPER) || defined(FASTTV)
        if (fastBeta != FLOAT_ZERO && ii == 0)
#if defined(FASTNLM)
            rhs += fastBeta * NLMGradient(
#ifdef FASTNLMLOCAL
                // Coordinates of the voxel within the local memory cache filled at the start
                lCacheF, d_gaussian, CINT(LID0) + SWINDOWX + PWINDOWX, CINT(LID1) + SWINDOWY + PWINDOWY, zz + SWINDOWZ + PWINDOWZ,
#else
                d_imNLM, d_gaussian, i.x, i.y, CINT(ind) + fastPriorZOffset,
#endif
                fastNLMh, fastEpps
#if NLTYPE >= 3
                , fastNLMgamma
#endif
#if NLTYPE == 6
                , fastNLMp, fastNLMq, fastNLMc
#endif
#ifdef NLMADAPTIVE
                , fastNLMs
#endif
#ifdef NLMREF
#ifdef FASTNLMLOCAL
                , lCacheRefF
#else
                , d_urefNLM
#endif
#endif
            );
#elif defined(FASTRDP)
#ifdef FASTRDPCORNERS
#ifdef FASTSWLOCAL
            rhs += fastBeta * priorGradientSW(lCacheSW, d_swWeight, CINT(LID0) + SWINDOWX, CINT(LID1) + SWINDOWY, zz + SWINDOWZ, fastEpps, fastRDPgamma);
#else
            rhs += fastBeta * priorGradientSW(d_imNLM, d_swWeight, i.x, i.y, CINT(ind) + fastPriorZOffset, fastEpps, fastRDPgamma
#ifdef PRIORREF
                , d_urefNLM
#endif
            );
#endif
#else
            rhs += fastBeta * RDPGradientNorm(d_imNLM, i.x, i.y, CINT(ind) + fastPriorZOffset, fastRDPgamma, fastEpps);
#endif
#elif defined(FASTGGMRF)
#ifdef FASTSWLOCAL
            rhs += fastBeta * priorGradientSW(lCacheSW, d_swWeight, CINT(LID0) + SWINDOWX, CINT(LID1) + SWINDOWY, zz + SWINDOWZ, fastEpps, fastGGMRFp, fastGGMRFq, fastGGMRFc, fastGGMRFpqc);
#else
            rhs += fastBeta * priorGradientSW(d_imNLM, d_swWeight, i.x, i.y, CINT(ind) + fastPriorZOffset, fastEpps, fastGGMRFp, fastGGMRFq, fastGGMRFc, fastGGMRFpqc);
#endif
#elif defined(FASTHYPER)
#ifdef FASTSWLOCAL
            rhs += fastBeta * priorGradientSW(lCacheSW, d_swWeight, CINT(LID0) + SWINDOWX, CINT(LID1) + SWINDOWY, zz + SWINDOWZ, fastEpps, fastHYPERsigma);
#else
            rhs += fastBeta * priorGradientSW(d_imNLM, d_swWeight, i.x, i.y, CINT(ind) + fastPriorZOffset, fastEpps, fastHYPERsigma);
#endif
#elif defined(FASTTV)
            // Standard/SATV/JPTV only
            rhs += fastBeta * TVGradient(d_imNLM, i.x, i.y, CINT(ind) + fastPriorZOffset, fastTVsigma, fastEpps);
#endif // FASTNLM / FASTRDP / FASTGGMRF / FASTHYPER / FASTTV dispatch
#endif // defined(FASTNLM) || defined(FASTRDP) || defined(FASTGGMRF) || defined(FASTHYPER) || defined(FASTTV)
        // Image-based preconditioning
#if FASTPRECOND == 1
        // Diagonal normalization preconditioner (type 0)
        rhs /= d_precond[idx];
#elif FASTPRECOND == 2
        // EM preconditioner (type 1)
#if defined(MBSREM) && FASTALG == 1
        // MBSREM feeds the EM preconditioner the image mirrored about U/2 (alpha == U here)
        const float precEst = (curEst >= fastAlpha * FLOAT_HALF) ? (fastAlpha - curEst) : curEst;
#else
        const float precEst = curEst;
#endif
        rhs *= precEst / d_precond[idx];
#endif
        // PDHG primal step (shared with auxKernels.cl PDHGUpdate), or the Poisson update (shared with
        // auxKernels.cl PoissonUpdate) for PKMA/MBSREM/BSREM
#if FASTALG == 0
        d_OSEM[idx + fastImOffset] = PDHGSubsetPrimal(curEst, rhs, fastTau, fastEpps, fastPositivity);
#else
        d_OSEM[idx + fastImOffset] = PoissonUpdateVoxel(curEst, rhs, fastLambda, fastEpps, fastAlpha, fastPositivity);
#endif
        if (no_norm == 0u)
            d_Summ[idx] = wSum[zz];
        idx += d_N.y * d_N.x;
    }
}

#endif // FASTPDHG
