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
    float ellipsePower = scalarParams.ellipsePower; \
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
    const float r = scalarParams.helicalRadius; \
    const float dL = scalarParams.dL; \
    const float global_factor = scalarParams.global_factor; \
    const float sigma_x = scalarParams.sigma_x; \
    const float d_epps = scalarParams.epps; \
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
    const int aa = scalarParams.aa; \
    const float bzGlobalMin = scalarParams.dSize5.x; \
    const float bzGlobalMax = scalarParams.dSize5.y;

#define UNPACK_SCALAR_PARAMS_4_BP(scalarParams) \
    const uint d_size_x = scalarParams.nRowsD; \
    const uint d_sizey = scalarParams.nColsD; \
    const float2 d_dPitch = scalarParams.dPitch; \
    const float r = scalarParams.helicalRadius; \
    const uint3 d_N = scalarParams.d_N; \
    const float3 b = scalarParams.b; \
    const float3 d_d = scalarParams.d; \
    const float2 dSize5 = scalarParams.dSize5; \
    const float kerroin = scalarParams.kerroin4; \
    const float DSC = scalarParams.DSC; \
    const uchar no_norm = scalarParams.no_norm; \
    const long d_nProjections = scalarParams.nProjections; \
    const int ii = scalarParams.aa; \
    const float bzGlobalMin = scalarParams.dSize5.x; \
    const float bzGlobalMax = scalarParams.dSize5.y;

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

#define UNPACK_SCALAR_PARAMS_5_BP(scalarParams) \
    const uint d_nRows = scalarParams.nRowsD; \
    const uint d_nCols = scalarParams.nColsD; \
    const float2 d_dPitch = scalarParams.dPitch; \
    const uint3 d_N = scalarParams.d_N; \
    const float3 b = scalarParams.b; \
    const float3 d_d = scalarParams.d; \
    const float3 d_scale = scalarParams.d_Scale5; \
    const float2 d_Size = scalarParams.dSize5; \
    const long d_nProjections = scalarParams.nProjections; \
    const uchar no_norm = scalarParams.no_norm; \
    const int ii = scalarParams.aa;

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

using metal::dot;

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
#define ATAN2 metal::atan2
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
#define CUINT_rte(a) static_cast<uint>(metal::rint((a)))
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
#define TEX10 [[texture(10)]]
#define TEX19 [[texture(19)]]
#else
#define TEX1 [[buffer(1)]]
#define TEX2 [[buffer(2)]]
#define TEX3 [[buffer(3)]]
#define TEX4 [[buffer(4)]]
#define TEX7 [[buffer(7)]]
#define TEX8 [[buffer(8)]]
#define TEX9 [[buffer(9)]]
#define TEX10 [[buffer(10)]]
#define TEX19 [[buffer(19)]]
#endif

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
// Fast math and the discrete model discontinuities (measured 2026-10-19): projector types 2 and 3 contain discontinuities in the model: whole
// voxel columns entering or leaving the ray for LORs grazing the corner of the image box (type 2), the ORTH slice TOF weights depending on the
// Siddon x/y crossing order at exact vertex ties (type 3) and the 0 -> THR step at the ORTH tube boundary (type 3). Fast division (CUDA
// --use_fast_math/__fdividef, OpenCL native_divide/-cl-fast-relaxed-math) and the implicit FMA contraction flip rays across these steps, so
// CUDA and OpenCL (and fast and strict builds) can differ by up to ~1e-2 of the maximum in a few elements (median ~5e-7). Precise division
// in the 6 Siddon pre-loop divisions and the 3 ORTH divisions together with no FMA contraction (FP_CONTRACT OFF, --fmad=false, no
// -cl-fast-relaxed-math) removes this at ~+3-7 % per FP+BP pair (types 2/3 FP +10-17 %). It is deliberately left as is.
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
// See the note on fast math and the discrete model discontinuities at the OpenCL DIVIDE definition
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
// Widened to 64-bit (via ULONG) to match OpenCL's get_global_id/get_global_size, which return
// size_t
// This keeps global-index arithmetic shared with OpenCL kernels (e.g. idx = GID0 + GID1 *
// d_size_x + GID2 * d_sizey * d_size_x) from silently overflowing 32-bit at large volumes
#define GID0 ((ULONG)threadIdx.x + (ULONG)blockIdx.x * (ULONG)blockDim.x)
#define GID1 ((ULONG)threadIdx.y + (ULONG)blockIdx.y * (ULONG)blockDim.y)
#define GID2 ((ULONG)threadIdx.z + (ULONG)blockIdx.z * (ULONG)blockDim.z)
#define GSIZE0 ((ULONG)blockDim.x * (ULONG)gridDim.x)
#define GSIZE1 ((ULONG)blockDim.y * (ULONG)gridDim.y)
#define GSIZE2 ((ULONG)blockDim.z * (ULONG)gridDim.z)
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

// Register budget of the TOF projector kernels, CUDA only (__launch_bounds__ means something else on AMD and is not used elsewhere):
// the multi-bin TOF kernels need many registers (up to ~170 with 33 bins) which can limit the occupancy to one block per multiprocessor.
// TOF_MAXREGS is the budget in registers per thread, the kernels ask for the number of resident blocks that fits into the
// 65536 registers of a multiprocessor. TOF_MAXREGS 0 disables the bound. Not defined (empty) without TOF
#if defined(CUDA) && !defined(HIP) && defined(TOF) && NBINS > 1
#ifndef TOF_MAXREGS
// Measured (RTX 4090, 256 threads per block, 5-33 bins): the kernels are fastest without a bound except for 23-29 bins, where the
// unbounded kernel needs 146-168 registers (1 block per multiprocessor) and 128 registers (2 blocks) is 8-25 % faster, even with a
// small amount of spilling. A bound of 85 registers (3 blocks) spills heavily with more than 15 bins
#if NBINS >= 23 && NBINS <= 29
#define TOF_MAXREGS 128
#else
#define TOF_MAXREGS 0
#endif
#endif
#if TOF_MAXREGS > 0
#ifdef LOCAL_SIZE2
#define TOF_THREADS (LOCAL_SIZE * LOCAL_SIZE2)
#else
#define TOF_THREADS (LOCAL_SIZE)
#endif
#define TOF_BLOCKS (((65536 / TOF_MAXREGS) / TOF_THREADS) > 1 ? ((65536 / TOF_MAXREGS) / TOF_THREADS) : 1)
#define TOF_BOUNDS __launch_bounds__(TOF_THREADS, TOF_BLOCKS)
#else
#define TOF_BOUNDS
#endif
#else
#define TOF_BOUNDS
#endif

// Number of elements of the per-bin arrays of the projectors (ax and TOFWeights). Listmode data (except for the sensitivity image) has
// only the bin of the event, so these arrays have a single element (always index 0) instead of NBINS elements indexed with the bin of
// the event: an array that is indexed with a run-time value is put to the local (stack) memory of the thread, i.e. every voxel would
// read and write it from there
#if defined(TOF) && defined(LISTMODE) && !defined(SENS)
#define TOF_NARR 1
#else
#define TOF_NARR NBINS
#endif

#ifdef TOF //////////////// TOF ////////////////
// Exact bin-integrated TOF weights
//
// For an emission at the TOF coordinate x (signed distance from the LOR center, measured in the same way as the bin centers
// TOFCenter) the probability of the event being recorded in the bin with the edges [eLo, eHi] is
//     p(x) = Phi((x - eLo) / sigma) - Phi((x - eHi) / sigma),
// where Phi is the standard normal CDF. TOFWeights[b] is the mean of p(x) over the ray segment (element) that starts at D and
// ends at D - element * DDsign, i.e. the exact (not sampled) bin probability averaged over the segment. The weights are NOT
// normalized over the bins: the sum over the bins is 1 when the TOF window covers the segment and smaller when it does not,
// since events outside the window are not recorded.
//
// With G(u) = u * Phi(u) + phi(u) (G' = Phi) the segment integral is
//     sigma * [G((x1 - eLo) / sigma) - G((x0 - eLo) / sigma) - G((x1 - eHi) / sigma) + G((x0 - eHi) / sigma)].
// For float precision G(u) is split into max(u, 0) + H(|u|), where H(t) = phi(t) - t * Q(t) (Q = upper tail) is the exponentially
// small tail part. The max(u, 0) terms of the four G values sum to the piecewise linear clamp expression (s1 - s0) in
// TOFEdge, which cancels exactly when both ends of the segment are on the same side of the bin. H is set to zero for
// t >= TOF_CUTOFF (it is below 2e-10 there), so the bins that are far from the segment get exactly zero weights. The H terms
// at the start of the segment are the H terms at the end of the previous segment, so they are cached (TOFF, see TOFInit and
// TOFEdge).
// The difference of the antiderivatives has an absolute error of about sigma * (float epsilon) whatever the length of the segment,
// i.e. the error of length * weight is about 1e-6 mm. The weights of very short segments are noisy (they are clamped to [0, 1]), but
// such segments contribute almost nothing, so there is no special treatment of them (except for the ORTH builds, see TOFSetup).
//
// The bin width is the (uniform) spacing of the bin centers, |TOFCenter[1] - TOFCenter[0]|. The centers are always in the
// center-first order (0, -w, +w, -2w, ...), with the possible TOF offset added to all of them, so the spacing is the same.
// For NBINS == 1 all events are accepted (weight 1)
//
// No erf/erfc is used anywhere (MSL has neither): H is computed with a polynomial fit (TOFH, TOFEdge below)
#define TOF_CUTOFF 6.f // Distance in sigmas after which the Gaussian is considered to be zero
#define TOF_MIDPOINT_LIMIT 1e-2f // ORTH builds only: segments shorter than this (in sigmas) use the midpoint value
#define TOF_INV_SQRT_2PI 0.3989422804f

// Clipped traversal of the listmode rays. A listmode event has a single TOF bin b (TOFIndex) and its weights are exactly zero when the
// whole voxel segment is farther than TOF_CUTOFF sigmas from the bin (see above), i.e. only the part of the ray inside the window
//     c_b - w / 2 - TOF_CUTOFF * sigma <= x <= c_b + w / 2 + TOF_CUTOFF * sigma,
// x the TOF coordinate (TOFDis) and w the bin width, can contribute. With TOF_CLIP the projectors do not visit the voxels outside of
// this window (the voxels are not cut), the results change only by floating point rounding. Everything that needs the whole ray
// (the length of the ray, the normalization of the backprojection and the attenuation) is still computed from the whole ray.
// Not used for the sensitivity image (all bins), the orthogonal and volume projectors (the weights are shared by the voxels of the
// tube around the ray), the single-bin case and CT/SPECT
#if defined(LISTMODE) && !defined(SENS) && NBINS > 1 && !defined(ORTH) && !defined(CT) && !defined(SPECT)
#define TOF_CLIP
#endif

DEVICE void TOFDis(const float tStart, const float LL, float* D, float* DD) {
	*D = LL * tStart - LL / 2.f;
	*DD = *D;
}

#ifdef TOF_CLIP
// The window [*xLo, *xHi] of the TOF coordinate of the bin TOFIndex of a listmode event outside of which the TOF weights are zero
DEVICE void TOFClipWindow(CONSTANT float* TOFCenter, const float invSigma, const int TOFIndex, PTR_THR float* xLo, PTR_THR float* xHi) {
	const float halfW = 0.5f * FABS(TOFCenter[1] - TOFCenter[0]);
	const float cut = TOF_CUTOFF / invSigma;
	*xLo = TOFCenter[TOFIndex] - halfW - cut;
	*xHi = TOFCenter[TOFIndex] + halfW + cut;
}
#endif

// H(t) = G(-t) = phi(t) - t * Q(t), t >= 0
// Computed as phi(t) * g(t), where g(t) = 1 - t * Q(t) / phi(t) = s^2 * (c0 + c1 s + ... + c6 s^6), s = 1 / (1 + 0.3 t), is a
// minimax-type fit on [0, 6]. The fit itself has a maximum absolute error of H of 8e-10 in double precision, but with the
// stored nine-digit coefficients it is 3.5e-9 (at t ~ 0.02, mostly from the rounding of the coefficient sum), and the float
// evaluation error of H is about 2e-7. This needs one exponential and one reciprocal instead of an exponential and a
// complementary error function
DEVICE float TOFH(const float t) {
	const float s = 1.f / (1.f + 0.3f * t);
	float g = FMAD(s, 0.278036759f, -1.48782828f);
	g = FMAD(s, g, 2.71388525f);
	g = FMAD(s, g, -1.70602154f);
	g = FMAD(s, g, 1.12047934f);
	g = FMAD(s, g, -0.029800149f);
	g = FMAD(s, g, 0.111248612f);
	// Exactly zero after TOF_CUTOFF sigmas (a select, not a branch)
	return ((t < TOF_CUTOFF) ? TOF_INV_SQRT_2PI * EXP(-0.5f * t * t) : 0.f) * s * s * g;
}

// The tail part of the antiderivative at the position x for a bin with the edges eLo and eHi
DEVICE float TOFTail(const float x, const float eLo, const float eHi, const float invSigma) {
	return TOFH(FABS((x - eLo) * invSigma)) - TOFH(FABS((x - eHi) * invSigma));
}

// Bin layout: the bins are on a uniform grid with the spacing dw = TOFCenter[1] - TOFCenter[0] (signed) in the center-first
// order produced by the host code, i.e. bin b has the center c0 + k_b * dw with k_b = (b + 1) / 2 for odd b and -b / 2 for even b
// (0, 1, -1, 2, -2, ...). Adjacent bins share an edge, so the NBINS bins have only NBINS + 1 unique edges
// P_e = c0 + (e + TOF_KMIN - 0.5) * dw, e = 0 ... NBINS, and the bin with k_b = e + TOF_KMIN lies between the edges e and e + 1.
// Everything (the H terms, the clamp terms) is evaluated once per edge. Since k_b is known at compile time, the arrays are
// indexed with compile-time constants after the loops are unrolled (no local memory); a layout that is not known at compile time
// would need dynamic indexing of per-thread arrays. TOFInit poisons the cache with NaN if TOFCenter does not follow this layout
// (the weights, and thus the output, become NaN), so a different layout cannot silently give wrong weights.
// For one edge the quantity a_e = (max(x1, P_e) - max(x0, P_e)) / sigma + H_e(x1) - H_e(x0) is formed, and the bin weight is
// sgn(dw) * (a_e - a_(e+1)) / ((x1 - x0) / sigma).
#define TOF_KMIN (-((NBINS - 1) / 2))
#define TOF_BIN(e) ((((e) + TOF_KMIN) > 0) ? (2 * ((e) + TOF_KMIN) - 1) : (-2 * ((e) + TOF_KMIN)))

#if NBINS > 1
#if defined(LISTMODE) && !defined(SENS)
#ifdef ORTH
// The quantity a_e of the edge at P for the segments shorter than TOF_MIDPOINT_LIMIT sigmas (the 'mid' case of TOFEdge, the same
// arithmetic): a_e = -Q(u_m), u_m = the middle of the segment in sigmas from the edge. halfdu is half the length of the segment in sigmas
DEVICE float TOFMidEdge(const float x1, const float P, const float invSigma, const float halfdu) {
	const float u = (x1 - P) * invSigma;
	const float t = FABS(u);
	const float s = 1.f / (1.f + 0.3f * t);
	const float E = (t < TOF_CUTOFF) ? TOF_INV_SQRT_2PI * EXP(-0.5f * t * t) : 0.f;
	float pq = FMAD(s, -0.155944471f, 0.384600595f);
	pq = FMAD(s, pq, 0.21560697f);
	pq = FMAD(s, pq, 0.153780407f);
	pq = FMAD(s, pq, 0.366935985f);
	pq = FMAD(s, pq, 0.288334552f);
	const float q = E * s * pq;
	return -(((u >= 0.f) ? q : 1.f - q) + E * halfdu);
}
#endif

// Listmode data (except for the sensitivity image) needs only the bin TOFIndex: both its edges are evaluated directly from
// TOFCenter[TOFIndex], hc is the cached H terms of the bin at x0 and is updated to x1
// The orthogonal and volume projectors (ORTH) use the same rule for short segments as the multi-bin path (TOFSetup): for segments
// shorter than TOF_MIDPOINT_LIMIT sigmas (including empty ones) the weight is the probability at the middle of the segment, the bin
// weight is then a(eLo) - a(eHi) (the same value as sgn(dw) * (a_e - a_(e+1)) of the multi-bin path for either sign of dw)
DEVICE float TOFSingle(const float x0, const float x1, CONSTANT float* TOFCenter, const int TOFIndex, const float invSigma, float* hc) {
	const float halfW = 0.5f * FABS(TOFCenter[1] - TOFCenter[0]);
	const float eLo = TOFCenter[TOFIndex] - halfW;
	const float eHi = TOFCenter[TOFIndex] + halfW;
	const float cut = TOF_CUTOFF / invSigma;
	float wt = 0.f;
	float F1 = 0.f;
	if (FMAX(x0, x1) >= eLo - cut && FMIN(x0, x1) <= eHi + cut) {
		F1 = TOFTail(x1, eLo, eHi, invSigma);
		const float du = (x1 - x0) * invSigma;
#ifdef ORTH
		if (FABS(du) < TOF_MIDPOINT_LIMIT) {
			const float halfdu = 0.5f * du;
			wt = TOFMidEdge(x1, eLo, invSigma, halfdu) - TOFMidEdge(x1, eHi, invSigma, halfdu);
			wt = (wt < 0.f) ? 0.f : ((wt > 1.f) ? 1.f : wt);
		}
		else
#endif
		if (du != 0.f) {
			const float s0 = FMIN(FMAX(x0, eLo), eHi);
			const float s1 = FMIN(FMAX(x1, eLo), eHi);
			wt = ((s1 - s0) * invSigma + (F1 - *hc)) / du;
			wt = (wt < 0.f) ? 0.f : ((wt > 1.f) ? 1.f : wt);
		}
	}
	*hc = F1;
	return wt;
}
#else
// Moves D to the end of the segment and computes the common values of the segment: the start x0, end x1, the first edge P0, the
// grid spacing dw and the scaling of the weights, sign(dw) / ((x1 - x0) / sigma) (zero for an empty segment: its weights are zero)
// The orthogonal and volume projectors (ORTH) apply the weights of a segment to all the voxels of the tube around it without
// scaling them with the length of the segment, so the weights of very short segments must be accurate as well: for segments shorter
// than TOF_MIDPOINT_LIMIT sigmas those builds use the probability at the middle of the segment (halfdu, mid)
DEVICE void TOFSetup(const float DDsign, const float element, CONSTANT float* TOFCenter, const float invSigma, float* D,
	float* x0, float* x1, float* P0, float* dw, float* scale
#ifdef ORTH
	, float* halfdu, int* mid
#endif
) {
	*x0 = *D;
	*x1 = FMAD(-element, DDsign, *x0);
	*D = *x1;
	const float c0 = TOFCenter[0];
	*dw = TOFCenter[1] - c0;
	*P0 = FMAD(*dw, -0.5f + TOF_KMIN, c0);
	const float du = (*x1 - *x0) * invSigma;
	const float sgnw = (*dw > 0.f) ? 1.f : -1.f;
#ifdef ORTH
	*halfdu = 0.5f * du;
	*mid = (FABS(du) < TOF_MIDPOINT_LIMIT) ? 1 : 0;
	*scale = (*mid) ? sgnw : sgnw / du;
#else
	*scale = (du != 0.f) ? sgnw / du : 0.f;
#endif
}

// The quantity a_e of the edge at P (see above). hc is the cached H(|x0 - P| / sigma), updated to the value at x1.
// There is no special case for short segments (except in the ORTH builds, see TOFSetup): the difference of the antiderivatives has an
// absolute error of about sigma * (float epsilon), independent of the segment length, so the error of length * weight is about
// 1e-6 mm; the weights of very short segments can be noisy (they are clamped to [0, 1]) but they contribute almost nothing
DEVICE float TOFEdge(const float x0, const float x1,
#ifdef ORTH
	const float halfdu, const int mid,
#endif
	const float P, const float invSigma, float* hc) {
	const float u = (x1 - P) * invSigma;
	const float t = FABS(u);
	const float s = 1.f / (1.f + 0.3f * t);
	float pg = FMAD(s, 0.278036759f, -1.48782828f);
	pg = FMAD(s, pg, 2.71388525f);
	pg = FMAD(s, pg, -1.70602154f);
	pg = FMAD(s, pg, 1.12047934f);
	pg = FMAD(s, pg, -0.029800149f);
	pg = FMAD(s, pg, 0.111248612f);
	// phi(t), exactly zero after TOF_CUTOFF sigmas (a select, not a branch): the bins that are far from the segment get exactly zero
	// weights and the H terms cancel exactly
	const float E = (t < TOF_CUTOFF) ? TOF_INV_SQRT_2PI * EXP(-0.5f * t * t) : 0.f;
	const float h1 = E * s * s * pg;
	float a;
#ifdef ORTH
	if (mid) {
		// a_e = -Q(u_m), u_m = the middle of the segment in sigmas from the edge, from the quantities of the end point with a first
		// order Taylor expansion (Q' = -phi, the second order term is below 3e-6 for these segments). Q(t) = phi(t) * s * (q0 + q1 s +
		// ... + q5 s^5) is a polynomial fit with a maximum absolute error of 4e-8 (double precision), sharing the exponential and s
		float pq = FMAD(s, -0.155944471f, 0.384600595f);
		pq = FMAD(s, pq, 0.21560697f);
		pq = FMAD(s, pq, 0.153780407f);
		pq = FMAD(s, pq, 0.366935985f);
		pq = FMAD(s, pq, 0.288334552f);
		const float q = E * s * pq;
		a = -(((u >= 0.f) ? q : 1.f - q) + E * halfdu);
	}
	else
#endif
		a = (FMAX(x1, P) - FMAX(x0, P)) * invSigma + (h1 - *hc);
	*hc = h1;
	return a;
}
#endif
#endif

// Initializes the cache of the H terms at the current position D (TOFF, NBINS + 1 values, or one value for the listmode bin
// TOFIndex). This has to be called every time D is set (the start of each ray) before the first call to the functions below
DEVICE void TOFInit(const float D, CONSTANT float* TOFCenter, const float invSigma, float* TOFF
#if defined(LISTMODE) && !defined(SENS)
	, const int TOFIndex
#endif
) {
#if NBINS > 1
#if defined(LISTMODE) && !defined(SENS)
	const float halfW = 0.5f * FABS(TOFCenter[1] - TOFCenter[0]);
	TOFF[0] = TOFTail(D, TOFCenter[TOFIndex] - halfW, TOFCenter[TOFIndex] + halfW, invSigma);
#else
	const float c0 = TOFCenter[0];
	const float dw = TOFCenter[1] - c0;
	// Check the layout of the bins
	float bad = 0.f;
#if !defined(__CUDACC__) && !defined(__HIPCC__)
#pragma unroll NBINS
#else
#pragma unroll
#endif
	for (int b = 1; b < NBINS; b++)
		bad += (FABS(TOFCenter[b] - FMAD(dw, CFLOAT((b & 1) ? ((b + 1) / 2) : (-b / 2)), c0)) > 1e-3f * FABS(dw)) ? 1.f : 0.f;
	// NaN if the layout is wrong, otherwise zero
	const float poison = SQRT(-bad);
	const float P0 = FMAD(dw, -0.5f + TOF_KMIN, c0);
#if !defined(__CUDACC__) && !defined(__HIPCC__)
#pragma unroll NBINS
#else
#pragma unroll
#endif
	for (int e = 0; e <= NBINS; e++)
		TOFF[e] = TOFH(FABS((D - FMAD(dw, CFLOAT(e), P0)) * invSigma)) + poison;
#endif
#endif
}

// Computes the weights of the bins for the segment starting at *D and ending at *D - element * DDsign, adjacent segments share
// an end point so the cache TOFF is kept up to date. This function also moves D to the end of the segment, i.e. D must not be
// modified elsewhere (after TOFInit). Listmode data (except for the sensitivity image) needs only the bin TOFIndex.
// Note that the weights are consumed immediately by denominator/rhs (and reused for several voxels by the orthogonal and volume
// projectors); the compiler keeps the unrolled arrays in registers, writing the weights and using them directly in the same loop
// was measured to give the same speed
DEVICE void TOFLoop(const float DDsign, const float element, CONSTANT float* TOFCenter, const float invSigma, float* D, float* TOFWeights, float* TOFF
#if defined(LISTMODE) && !defined(SENS)
	, const int TOFIndex
#endif
) {
#if NBINS == 1
	*D = FMAD(-element, DDsign, *D);
	TOFWeights[0] = 1.f;
#elif defined(LISTMODE) && !defined(SENS)
	const float x0 = *D;
	const float x1 = FMAD(-element, DDsign, x0);
	*D = x1;
	float hc = TOFF[0];
	TOFWeights[0] = TOFSingle(x0, x1, TOFCenter, TOFIndex, invSigma, &hc);
	TOFF[0] = hc;
#else
	float x0, x1, P0, dw, scale;
#ifdef ORTH
	float halfdu;
	int mid;
	TOFSetup(DDsign, element, TOFCenter, invSigma, D, &x0, &x1, &P0, &dw, &scale, &halfdu, &mid);
	float aPrev = TOFEdge(x0, x1, halfdu, mid, P0, invSigma, &TOFF[0]);
#else
	TOFSetup(DDsign, element, TOFCenter, invSigma, D, &x0, &x1, &P0, &dw, &scale);
	float aPrev = TOFEdge(x0, x1, P0, invSigma, &TOFF[0]);
#endif
#if !defined(__CUDACC__) && !defined(__HIPCC__)
#pragma unroll NBINS
#else
#pragma unroll
#endif
	for (int e = 0; e < NBINS; e++) {
#ifdef ORTH
		const float aNext = TOFEdge(x0, x1, halfdu, mid, FMAD(dw, CFLOAT(e + 1), P0), invSigma, &TOFF[e + 1]);
#else
		const float aNext = TOFEdge(x0, x1, FMAD(dw, CFLOAT(e + 1), P0), invSigma, &TOFF[e + 1]);
#endif
		// A probability: clamped to [0, 1]. Not FMAX/FMIN: a NaN (wrong bin layout, see TOFInit) must not be turned into a number
		const float wt = scale * (aPrev - aNext);
		TOFWeights[TOF_BIN(e)] = (wt < 0.f) ? 0.f : ((wt > 1.f) ? 1.f : wt);
		aPrev = aNext;
	}
#endif
}
#endif //////////////// END TOF ////////////////

#if defined(N_RAYS)
DEVICE void multirayCoordinateShiftXY(PTR_THR FLOAT3 *s, PTR_THR FLOAT3 *d, const int lor, const float cr) {
	const float dx = (*d).x - (*s).x;
	const float dy = (*d).y - (*s).y;
	const float lenXY = SQRT(dx * dx + dy * dy);
	if (lenXY < 1e-8f)
		return;
	const float interval = ((CFLOAT(lor) + FLOAT_HALF) * (cr / CFLOAT(N_RAYS2D)) - cr * FLOAT_HALF) / lenXY;
	const float shiftX = -dy * interval;
	const float shiftY = dx * interval;
	(*s).x += shiftX;
	(*d).x += shiftX;
	(*s).y += shiftY;
	(*d).y += shiftY;
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
	, float* TOFWeights
#endif
) {
	float apu = 0.f;
	forwardProject(local_ele, &apu, localInd, d_OSEM);
#ifdef TOF
#if defined(LISTMODE) && !defined(SENS)
	const int to = 0;
#else
#if !defined(__CUDACC__) && !defined(__HIPCC__)
#pragma unroll NBINS
#endif
	for (int to = 0; to < NBINS; to++) {
#endif
		ax[to] += TOFWeights[to] * apu;
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
	, float* TOFWeights
#endif
) {
#ifdef TOF
	FLOAT val = FLOAT_ZERO;
	FLOAT yaxTOF = FLOAT_ZERO;
#if defined(LISTMODE) && !defined(SENS)
	const int to = 0;
#else
#if !defined(__CUDACC__) && !defined(__HIPCC__)
#pragma unroll NBINS
#endif
	for (int to = 0; to < NBINS; to++) {
#endif
		const float apu = local_ele * TOFWeights[to];
		val += apu;
		yaxTOF += (apu * ax[to]);
#if !defined(LISTMODE) || defined(SENS)
	}
#endif
#else
	float yaxTOF = ax[0] * local_ele;
	float val = local_ele;
#endif
	// test
	if (ISNAN(val))
		val = 0.f;
	if (ISNAN(yaxTOF))
		yaxTOF = 0.f;
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
#if defined(N_RAYS) && !defined(PROJ5)
	, const int lorXY, const int lorZ
#endif
#ifdef PROJ5
	, PTR_THR float3* dR, PTR_THR float3* dL, PTR_THR float3* dU, PTR_THR float3* dD
#endif
) {
	int id = i.z * 6;
	*s = CMFLOAT3(d_xyz[id], d_xyz[id + 1], d_xyz[id + 2]);
	*d = CMFLOAT3(d_xyz[id + 3], d_xyz[id + 4], d_xyz[id + 5]);
	float2 indeksi = MFLOAT2(CFLOAT(i.x) - CFLOAT(d_size_x) / 2.f + .5f, CFLOAT(i.y) - CFLOAT(d_sizey) / 2.f + .5f);
#if defined(N_RAYS) && !defined(PROJ5)
	if (N_RAYS2D > 1)
		indeksi.x += (CFLOAT(lorXY) + .5f) / CFLOAT(N_RAYS2D) - .5f;
	if (N_RAYS3D > 1)
		indeksi.y += (CFLOAT(lorZ) + .5f) / CFLOAT(N_RAYS3D) - .5f;
#endif
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
    const float ellipsePower
) {
    // Box support. Hosts pass FLT_MAX instead of Inf since isinf() is unreliable
    // under fast-math (-ffast-math, -cl-fast-relaxed-math, Metal fast math)
    if (ellipsePower > 1.0e20f) {
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
    } else if (ellipsePower == 2.f) {
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
    const float ellipsePower
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
	// With multiple layers z has three values per sinogram, the third one is the layer combination of the sinogram. The
	// transaxial coordinates are stored as one block (of d_size_x * d_sizey bins) per layer combination while d_xyindex
	// is local to the block
#if defined(NLAYERS)
	const uint indz = d_zindex[idx] * 3;
	const uint layer = CUINT(d_z[indz + 2]);
	const uint ind = d_xyindex[idx] * 4 + layer * d_size_x * d_sizey * 4;
#else
	const uint indz = d_zindex[idx] * 2;
	const uint ind = d_xyindex[idx] * 4;
#endif
	(*s).x = d_xy[ind];
	(*s).y = d_xy[ind + 1];
	(*d).x = d_xy[ind + 2];
	(*d).y = d_xy[ind + 3];
	(*s).z = d_z[indz];
	(*d).z = d_z[indz + 1];
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
	, const uint d_sizey
#endif
) {
	const int id = (i.x + i.y * d_size_x) * 4;
#if defined(NLAYERS)
	// With multiple layers z has three values per sinogram, the third one is the layer combination of the sinogram
	// (the transaxial coordinates are stored as one block of d_size_x * d_sizey bins per layer combination)
	const int idz = i.z * 3;
	const int layer = CINT(d_z[idz + 2]);
	*s = CMFLOAT3(d_xy[id + layer * d_size_x * d_sizey * 4], d_xy[id + layer * d_size_x * d_sizey * 4 + 1], d_z[idz]);
	*d = CMFLOAT3(d_xy[id + layer * d_size_x * d_sizey * 4 + 2], d_xy[id + layer * d_size_x * d_sizey * 4 + 3], d_z[idz + 1]);
#else
	const int idz = i.z * 2;
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
// A coordinate that is constant along the ray (axis-parallel rays) is on a voxel boundary if it is closer than this to
// one (in voxel units). Which of the two neighbouring voxels floor() selects for such a ray is decided by rounding (it differs
// between the backends and the fast math modes), the upper voxel is thus always used. The error of the coordinate in voxel units is
// orders of magnitude smaller than this.
#define BOUNDARY_TOL 1e-4f

// Voxel index of a coordinate that is constant along the ray, u = (coordinate - b) / voxel size is the coordinate in voxel units:
// floor(u), except that a ray on a boundary between voxels k - 1 and k is in the voxel k. Not used by the orthogonal and volume-based
// projectors (ORTH), they use floor(u)
DEVICE int tieFloor(const float u) {
#ifdef ORTH
	return CINT(FLOOR(u));
#else
	const float r = FLOOR(u + 0.5f);
	return CINT((FABS(u - r) < BOUNDARY_TOL) ? r : FLOOR(u));
#endif
}

// Compute the voxel index where the current perpendicular measurement starts
DEVICE int perpendicular_start(const float d_b, const float d, const float d_d, const uint d_N) {
	int tempi = tieFloor((d - d_b) / d_d);
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
	const float local_norm, const float L, const float pathLength, const int apuX1, const int apuX2, const float dT1, const float dT2) {
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
#if defined(ORTH)
	float temp = FLOAT_ONE;
#elif defined(N_RAYS) //////////////// MULTIRAY ////////////////
#ifdef TOTLENGTH
	float temp = FLOAT_ONE / (L * CFLOAT(N_RAYS));
#else
	float temp = FLOAT_ONE / (pathLength * CFLOAT(N_RAYS));
#endif
#else
#ifdef TOTLENGTH
	float temp = FLOAT_ONE / L;
#else
	float temp = FLOAT_ONE / pathLength;
#endif
#endif //////////////// END MULTIRAY ////////////////
#if defined(ATN) && !defined(SPECT) //////////////// ATTENUATION ////////////////
		float jelppi = FLOAT_ZERO;
#ifdef USEIMAGES
		int3 atnind = *tempi;
#else
		LONG atnind = *z_loop;
#endif
		for (int iii = apuX1; iii <= apuX2; iii++) {
			float d_in = d_d2;
			if (iii == apuX1 && iii == apuX2)
				d_in = FMAX(dT1 + dT2 - d_d2, FLOAT_ZERO);
			else if (iii == apuX1)
				d_in = dT1;
			else if (iii == apuX2)
				d_in = dT2;
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
			compute_attenuation(d_in, atnind, d_atten, &jelppi, ii);
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
		const uint v = CUINT_rtp((p_t - b) / d);
		*v_min = v < 1u ? 1u : v;
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
	*t_0 = (b + CFLOAT(*v_min) * d - s) / diff;
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
		const uint v = CUINT_sat_rtz((p_t - b) / d);
		*v_max = v > N - 1u ? N - 1u : v;
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
	*t_0 = (b + CFLOAT(*v_max) * d - s) / diff;
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
#ifdef ORTH
	// The voxel [tc, t0] overlaps the line segment t = [0, 1] iff t0 >= 0 and tc <= 1 (this includes a segment that is
	// completely inside the voxel, where neither t0 nor tc is inside [0, 1])
	*pass = (*t0 >= 0.f && *tc <= 1.f);
	float local_ele = 0.f;
	if (*pass)
		local_ele = voxelValue(FMIN(*t0, 1.f), FMAX(*tc, 0.f), L);
#else
	// The part of the voxel segment [tc, t0] that is inside the line segment t = [0, 1], zero if the voxel is outside of it.
	// The overlap length is clipped to the line segment [0, 1], so a segment that is completely inside one voxel gets its
	// actual length. This is evaluated for every voxel of every ray, so it is written with min/max/selects only, without
	// the chain of comparisons that decides whether the voxel overlaps the segment
	const float local_ele = FMAX(FMIN(*t0, 1.f) - FMAX(*tc, 0.f), 0.f) * L;
	*pass = local_ele > 0.f;
#endif
	*temp_ijk += u;
	*tc = *t0;
	*t0 += tu;
	return local_ele;
}

// initial voxel index (beginning of the ray) from the same plane indices as the first crossing times,
// so both stay consistent under rounding
DEVICE int start_voxel(const int u, const uint v_min, const uint v_max, const uint N) {
	const int v = (u > 0) ? (int)v_min - 1 : (int)v_max;
	return v < 0 ? 0 : (v > (int)N - 1 ? (int)N - 1 : v);
}

DEVICE bool siddon_pre_loop_2D(const float b1, const float b2, const float diff1, const float diff2, const float max1, const float max2,
	const float d1, const float d2, const uint N1, const uint N2, PTR_THR int *temp1, PTR_THR int *temp2, PTR_THR float *t1u, PTR_THR float *t2u, PTR_THR uint *Np,
	const int TYYPPI, const float ys, const float xs, const float yd, const float xd, PTR_THR float *tc, PTR_THR int *u1, PTR_THR int *u2, PTR_THR float *t10, PTR_THR float *t20, PTR_THR bool *xy, PTR_THR float *tMaxOut) {
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
	*tMaxOut = tmax;
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

	// (26)
	*temp1 = start_voxel(*u1, imin, imax, N1);
	// (27)
	*temp2 = start_voxel(*u2, jmin, jmax, N2);

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
	PTR_THR float *ty0, PTR_THR float *tz0, PTR_THR bool *xy, const int3 ii, PTR_THR float *tMaxOut) {

	const float3 apuT = (float3)b - (float3)s;
	const float3 t0 = apuT / (float3)diff;
	const float3 tBack = DIVIDE3((float3)max - (float3)s, (float3)diff);

	const float3 tMin = FMIN(t0, tBack);
	const float3 tMax = FMAX(t0, tBack);

	*tc = FMAX(FMAX(tMin.x, tMin.z), tMin.y);
	const float tmax = FMIN(FMIN(tMax.x, tMax.z), tMax.y);
	*tMaxOut = tmax;
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

	*tempi = start_voxel(*i, imin, imax, N.x);
	*tempj = start_voxel(*j, jmin, jmax, N.y);
	*tempk = start_voxel(*k, kmin, kmax, N.z);

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
    const float3 s,
    const float3 v,
    const float xc, const float yc, const float r,
    const float cosA, const float sinA,
    const float zMin, const float zMax,
    const float thetaSpan,
	float* theta,
    float3* intersection
) {
	const float fx = s.x - xc;
	const float fy = s.y - yc;
    const float a = v.x * v.x + v.y * v.y;
    const float b = 2.f * (v.x * fx + v.y * fy);
    const float c = fx * fx + fy * fy - r * r;
    const float disc = b * b - 4.f * a * c;
    if (disc < 0.f || a == 0.f)
		return 0;

    for (int signv = -1; signv <= 1; signv += 2) {
        const float t = (-b + (float)signv * SQRT(disc)) / (2.f * a);
        if (t <= 0.f)
			continue;
		const float3 xyz = s + v * t;
		const float wx = xyz.x - xc;
		const float wy = xyz.y - yc;
		// Signed arc angle from the detector center, positive towards increasing row index (see getDetectorCoordinatesCT)
		const float th = ATAN2(-wx * sinA - wy * cosA, wx * cosA - wy * sinA);
		if (FABS(th) > thetaSpan / 2.f || xyz.z < zMin || xyz.z > zMax)
			continue;

        *intersection = xyz;
		*theta = th + thetaSpan / 2.f;
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
#ifdef NLMAXWEIGHT
    float maxWeight = FLOAT_ZERO;
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
#ifdef NLMAXWEIGHT
                maxWeight = FMAX(maxWeight, weight);
#endif
                const float uk = NLMFETCH(uIn, x + i, y + j, z + k);
                // Different NLM regularization methods
                // NLTYPE 0 = MRF NLM
                // NLTYPE 1 = NLTV
                // NLTYPE 2 = NLM filtered (i.e. similar to MRP)
                // NLTYPE 3 = NLRD
                // NLTYPE 4 = NL Lange
                // NLTYPE 5 = NLM filtered with Lange
                // NLTYPE 6 = NLGGMRF
                // NLTYPE 7 = NLGM
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
				// NLGM (Non-local Geman-McClure)
                const float delta = uj - uk;
                const float apu = (delta * delta + gamma * gamma);
                output += weight * (FLOAT_TWO * gamma * gamma * delta) / (apu * apu);
#else
                //NLTV
                const float apuU = uj - uk;
                output += (weight * apuU);
                outputAla += weight * apuU * apuU;
#endif // END NLM NLTYPE
            }
        }
    }
#ifdef NLMAXWEIGHT
#if NLTYPE == 2 || NLTYPE == 5
    output += maxWeight * uj;
#endif
    weight_sum += maxWeight;
#endif
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
// Buffer read with each coordinate clamped to the volume (clamp-to-edge, i.e. Neumann boundary), identical to the
// texture sampler used in the image version
DEVICE float SWBufRead(const CLGLOBAL float* CLRESTRICT buf, const int x, const int y, const int z, const int3 N) {
    const int xc = MAX(0, MIN(x, N.x - 1));
    const int yc = MAX(0, MIN(y, N.y - 1));
    const int zc = MAX(0, MIN(z, N.z - 1));
    return buf[xc + yc * N.x + zc * N.x * N.y];
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

// Buffer reads (SWBufRead) clamp to the volume like the texture sampler used in the image version

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

#ifdef ANATOMICAL3
// APLS vector q = (g - (xi.g) xi) / sqrt(|g|^2 - (xi.g)^2 + epps), xi = r / sqrt(|r|^2 + eta^2)
// g = forward differences of the image, r = forward differences of the reference image
DEVICE float3 APLSq(const float3 g, const float3 r, const float eta, const float epps) {
    const float rNorm = SQRT(r.x * r.x + r.y * r.y + r.z * r.z + eta * eta);
    const float xix = r.x / rNorm;
    const float xiy = r.y / rNorm;
    const float xiz = r.z / rNorm;
    const float gxi = g.x * xix + g.y * xiy + g.z * xiz;
    const float p = SQRT(FMAX(g.x * g.x + g.y * g.y + g.z * g.z - gxi * gxi, 0.f) + epps);
    return CMFLOAT3((g.x - gxi * xix) / p, (g.y - gxi * xiy) / p, (g.z - gxi * xiz) / p);
}
#endif

#ifdef ANATOMICAL1
// Anatomical TV type 1: q = 0.5 * grad_g Q(g) / sqrt(Q(g) + epps), Q(g) = g^T A g with the 3x3 weight matrix A of the voxel
// The 9 weights S[n + NN * kk] are read at the (clamped) voxel (x, y, z), g = forward differences of the image at that voxel
DEVICE float3 TVAnatQ(const CLGLOBAL float* CLRESTRICT S, const int x, const int y, const int z, const int3 N, const float3 g, const float epps) {
    const int xc = MAX(0, MIN(x, N.x - 1));
    const int yc = MAX(0, MIN(y, N.y - 1));
    const int zc = MAX(0, MIN(z, N.z - 1));
    const int NN = N.x * N.y * N.z;
    const int nc = xc + yc * N.x + zc * N.x * N.y;
    const float s0 = S[nc];
    const float s1 = S[nc + NN];
    const float s2 = S[nc + NN * 2];
    const float s3 = S[nc + NN * 3];
    const float s4 = S[nc + NN * 4];
    const float s5 = S[nc + NN * 5];
    const float s6 = S[nc + NN * 6];
    const float s7 = S[nc + NN * 7];
    const float s8 = S[nc + NN * 8];
    const float sxy = s1 + s3;
    const float sxz = s2 + s6;
    const float syz = s5 + s7;
    const float p = SQRT(s0 * g.x * g.x + s4 * g.y * g.y + s8 * g.z * g.z + sxy * g.x * g.y + sxz * g.x * g.z + syz * g.y * g.z + epps);
    return CMFLOAT3(0.5f * (2.f * s0 * g.x + sxy * g.y + sxz * g.z) / p,
        0.5f * (2.f * s4 * g.y + sxy * g.x + syz * g.z) / p,
        0.5f * (2.f * s8 * g.z + sxz * g.x + syz * g.y) / p);
}
#endif

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
    // phi_n = sqrt(Q_n(g_n) + epps) with the weights of voxel n, gradient = -(q_m.x + q_m.y + q_m.z) + q_{m-ex}.x + q_{m-ey}.y + q_{m-ez}.z
    // The image is clamped at the volume edge, so the differences that cross the edge are zero and do not contribute
    const float3 q0 = TVAnatQ(S, x, y, z, tvN, uijkP - uijk, epps);
    const float3 q1 = TVAnatQ(S, x - 1, y, z, tvN, u1, epps);
    const float3 q2 = TVAnatQ(S, x, y - 1, z, tvN, u2, epps);
    const float3 q3 = TVAnatQ(S, x, y, z - 1, tvN, u3, epps);
    float res = 0.f;
    if (x < tvN.x - 1)
        res -= q0.x;
    if (y < tvN.y - 1)
        res -= q0.y;
    if (z < tvN.z - 1)
        res -= q0.z;
    if (x > 0)
        res += q1.x;
    if (y > 0)
        res += q2.y;
    if (z > 0)
        res += q3.z;
    return res + 1e-7f;
#elif defined(ANATOMICAL2) // TV type 2
    // The reference image is clamped at the volume edge (same as the image)
    const float3 uijkR = MFLOAT3(SWBufRead(S, x + 1, y, z, tvN), SWBufRead(S, x, y + 1, z, tvN), SWBufRead(S, x, y, z + 1, tvN));
    const float3 apuS = (uijkR - SWBufRead(S, x, y, z, tvN));
    const float3 apu = uijkP - uijk;
    const float pvalijk = SQRT(apu.x * apu.x + apu.y * apu.y + apu.z * apu.z + C * (apuS.x * apuS.x + apuS.y * apuS.y + apuS.z * apuS.z) + epps);
    // Reference image differences at the three backward neighbors, needed for their phi = sqrt(|g|^2 + C |r|^2 + epps)
    const float vC = SWBufRead(S, x, y, z, tvN);
    const float vX = SWBufRead(S, x - 1, y, z, tvN);
    const float3 r1 = MFLOAT3(vC - vX, SWBufRead(S, x - 1, y + 1, z, tvN) - vX, SWBufRead(S, x - 1, y, z + 1, tvN) - vX);
    const float vY = SWBufRead(S, x, y - 1, z, tvN);
    const float3 r2 = MFLOAT3(SWBufRead(S, x + 1, y - 1, z, tvN) - vY, vC - vY, SWBufRead(S, x, y - 1, z + 1, tvN) - vY);
    const float vZ = SWBufRead(S, x, y, z - 1, tvN);
    const float3 r3 = MFLOAT3(SWBufRead(S, x + 1, y, z - 1, tvN) - vZ, SWBufRead(S, x, y + 1, z - 1, tvN) - vZ, vC - vZ);
    const float p1 = SQRT(u1.x * u1.x + u1.y * u1.y + u1.z * u1.z + C * (r1.x * r1.x + r1.y * r1.y + r1.z * r1.z) + epps);
    const float p2 = SQRT(u2.x * u2.x + u2.y * u2.y + u2.z * u2.z + C * (r2.x * r2.x + r2.y * r2.y + r2.z * r2.z) + epps);
    const float p3 = SQRT(u3.x * u3.x + u3.y * u3.y + u3.z * u3.z + C * (r3.x * r3.x + r3.y * r3.y + r3.z * r3.z) + epps);
    return (3.f * uijk - uijkP.x - uijkP.y - uijkP.z) / pvalijk + (uijk - uijkM.x) / p1 + (uijk - uijkM.y) / p2 + (uijk - uijkM.z) / p3 + 1e-7f;
#elif defined(ANATOMICAL3) // APLS
    // Asymmetric parallel level sets (Ehrhardt et al. 2016): P(u) = sum sqrt(|g|^2 - (xi.g)^2 + epps), xi = r / sqrt(|r|^2 + eta^2)
    // The gradient is the sum of the unit-like vectors q of the voxel itself and of its three backward neighbors
    // The reference image is clamped at the volume edge, so its differences vanish across the boundary
    const float vC = SWBufRead(S, x, y, z, tvN);
    const float3 r0 = MFLOAT3(SWBufRead(S, x + 1, y, z, tvN) - vC, SWBufRead(S, x, y + 1, z, tvN) - vC, SWBufRead(S, x, y, z + 1, tvN) - vC);
    const float vX = SWBufRead(S, x - 1, y, z, tvN);
    const float3 r1 = MFLOAT3(vC - vX, SWBufRead(S, x - 1, y + 1, z, tvN) - vX, SWBufRead(S, x - 1, y, z + 1, tvN) - vX);
    const float vY = SWBufRead(S, x, y - 1, z, tvN);
    const float3 r2 = MFLOAT3(SWBufRead(S, x + 1, y - 1, z, tvN) - vY, vC - vY, SWBufRead(S, x, y - 1, z + 1, tvN) - vY);
    const float vZ = SWBufRead(S, x, y, z - 1, tvN);
    const float3 r3 = MFLOAT3(SWBufRead(S, x + 1, y, z - 1, tvN) - vZ, SWBufRead(S, x, y + 1, z - 1, tvN) - vZ, vC - vZ);
    const float3 q0 = APLSq(uijkP - uijk, r0, C, epps);
    const float3 q1 = APLSq(u1, r1, C, epps);
    const float3 q2 = APLSq(u2, r2, C, epps);
    const float3 q3 = APLSq(u3, r3, C, epps);
    return -(q0.x + q0.y + q0.z) + q1.x + q2.y + q3.z + 1e-7f;
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
