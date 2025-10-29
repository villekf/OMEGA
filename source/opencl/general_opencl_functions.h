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
* Copyright (C) 2019-2024 Ville-Veikko Wettenhovi, Niilo Saarlemo
*
* This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************************************************************************/

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
#if defined(PTYPE4)
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

#ifdef METAL
#ifdef HALF // 16-bit floating point
#define CAST float
#define CFLOAT(a) static_cast<half>(a)
#define CFLOAT3(a) half3(a)
#define CMFLOAT3(x,y,z) half3((x), (y), (z))

#define make_float2(a,b) half2((a),(b))
#define make_float3(a,b,c) half3((a),(b),(c)) // TODO: replace with half
#define MFLOAT2(a,b) half2((a), (b)) // TODO: replace with half
#else // 32-bit floating point
#define CAST float
#define CFLOAT(a) static_cast<float>(a)
#define CFLOAT3(a) float3(a)
#define CMFLOAT3(x,y,z) float3((x), (y), (z))
#define FLOAT float
#define make_float2(a,b) float2((a),(b))
#define make_float3(a,b,c) float3((a),(b),(c))
#define MFLOAT2(a,b) float2((a), (b))
#endif

#define ACOS metal::acos
#define ALL metal::all
#define ANY metal::any
#define CINT(a)   static_cast<int>(a)
#define CINT_rtz(a) static_cast<int>(metal::trunc((a)))
#define CINT3_rtz(a) static_cast<int3>(a)
#define CLAMP3(a,b,c) metal::clamp((a),(b),(c))
#define CLGLOBAL device
#define CLONG_rtz(a) static_cast<long>(a)
#define CLRESTRICT
#define CMINT3(a,b,c) int3((a),(b),(c))
#define CONSTANT constant
#define CUINT(a) (uint)(a)
#define CUINT_rtp(a) static_cast<uint>(metal::ceil((a)))
#define CUINT_sat_rtz(a) static_cast<uint>(metal::clamp(metal::trunc(((float)a)), 0.0f, 4294967295.0f)) // TODO replace float with FLOAT
#define DEVICE inline
#define DIVIDE(a,b) ((a) / (b))
#define DIVIDE3(a,b) ((a) / (b))
#define EXP(a) metal::exp(a)
#define FABS metal::fabs
#define FMAD(a,b,c) metal::fma(a,b,c)
#define FMAD3(a,b,c) metal::fma((a),(b),(c))
#define FMAX metal::fmax
#define FMIN metal::fmin
#define ISINF metal::isinf
#define ISNAN metal::isnan
#define KERNEL kernel
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
#define MINT3(a,b,c) int3((a),(b),(c))
#define MUINT3(a, b, c) uint3(a, b, c)
#define POWR metal::pow
#define PTR_DEV device
#define PTR_THR thread
#define PTR_CONST constant
#define PTR_TG threadgroup
#define RCP(x) (1.f / x)
#define SQRT metal::sqrt

// Metal function definitions
inline FLOAT dot(half3 a, half3 b) {
    return metal::dot(a, b);
}
inline FLOAT dot(float3 a, float3 b) {
    return metal::dot(a, b);
}

inline void atomicAdd(volatile device metal::atomic_float* addr, float val)
{
    atomic_fetch_add_explicit(addr, val, metal::memory_order_relaxed);
}

// Metal scalar params unpacking
#define UNPACK_METAL_PARAMS(staticParams, dynamicParams) \
    FLOAT global_factor = staticParams.global_factor; \
    FLOAT d_epps = staticParams.epps; \
	uint d_size_x = staticParams.nRowsD; \
	uint d_det_per_ring = staticParams.det_per_ring; \
	FLOAT sigma_x = staticParams.sigma_x; \
	FLOAT coneOfResponseStdCoeffA = staticParams.coneOfResponseStdCoeffA; \
    FLOAT coneOfResponseStdCoeffB = staticParams.coneOfResponseStdCoeffB; \
    FLOAT coneOfResponseStdCoeffC = staticParams.coneOfResponseStdCoeffC; \
	FLOAT crystalSizeX = staticParams.dPitchX; \
	FLOAT crystalSizeY = staticParams.dPitchY; \
	FLOAT bmin = staticParams.bmin; \
	FLOAT bmax = staticParams.bmax; \
	FLOAT Vmax = staticParams.Vmax; \
	uint d_sizey = staticParams.nColsD; \
	long d_nProjections = dynamicParams.nProjections; \
	uint rings = staticParams.rings; \
	uint d_Nx = dynamicParams.d_N[0]; \
	uint d_Ny = dynamicParams.d_N[1]; \
	uint d_Nz = dynamicParams.d_N[2]; \
	FLOAT d_dx = dynamicParams.d[0]; \
	FLOAT d_dy = dynamicParams.d[1]; \
	FLOAT d_dz = dynamicParams.d[2]; \
	FLOAT bx = dynamicParams.b[0]; \
	FLOAT by = dynamicParams.b[1]; \
	FLOAT bz = dynamicParams.b[2]; \
	FLOAT d_bmaxx = dynamicParams.bmax[0]; \
	FLOAT d_bmaxy = dynamicParams.bmax[1]; \
	FLOAT d_bmaxz = dynamicParams.bmax[2]; \
	unsigned char no_norm = dynamicParams.no_norm; \
	unsigned long m_size = dynamicParams.m_size; \
	uint currentSubset = dynamicParams.currentSubset; \
	int aa = dynamicParams.aa; \
    float orthWidth = dynamicParams.orthWidth;
	
#endif
#ifdef OPENCL
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
#define POWR native_powr
#define POWN pown
#define RCP(x) native_recip(x)
#define FLOOR floor
#define CEIL ceil
#define ATAN2 atan2
#define ACOS acos
#define LOG native_log
#define CLAMP3(a, b, c) clamp(a, b, c)
#define CINT(a) convert_int(a)
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

#if defined(MASKBP) && defined(PTYPE4) && !defined(CT) && defined(BP)
__constant sampler_t sampler_MASK4 = CLK_NORMALIZED_COORDS_TRUE | CLK_FILTER_NEAREST | CLK_ADDRESS_CLAMP_TO_EDGE;
#else
__constant sampler_t sampler_MASK = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_CLAMP_TO_EDGE;
#endif
#endif
#if defined(CUDA)
#define PTR_DEV // Metal requires address space qualifier for pointers
#define PTR_THR 
#define PTR_CONST
#define PTR_TG
#define MIN min
#define FABS fabs
#define LENGTH length
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
#define CUINT_rtp(a) __float2uint_ru(a)
#define CUINT_rtz(a) __float2uint_rd(a)
#define CUINT_rte(a) __float2uint_rn(a)
#define CUINT_sat_rtz(a) __float2uint_rd(a)
#define CINT(a) (int)(a)
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
#define IMAGE3D cudaTextureObject_t
#define IMAGE2D cudaTextureObject_t
#define MUINT2(a, b) make_uint2(a, b)
#define MINT3(a, b, c) make_int3(a, b, c)
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
#define LOG logf

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

inline __device__ int3 operator-(int3 a, int3 b) {
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

#if (defined(MASKBP) && !defined(PTYPE4)) // This is due to projector type 4 using sampler_MASK4 in BP mask (but only in forward projection)
DEVICE int readMaskBP(MASKBPTYPE maskBP, typeT ind) {
    return 
    #ifdef USEIMAGES
    #ifdef CUDA
    #ifdef MASKBP3D
        tex3D<unsigned char>(maskBP, ind.x, ind.y, ind.z);
    #else
        tex2D<unsigned char>(maskBP, ind.x, ind.y);
    #endif
    #else
    #ifdef MASKBP3D
        read_imageui(maskBP, sampler_MASK, (int4)(ind.x, ind.y, ind.z, 0)).w;
    #else
        read_imageui(maskBP, sampler_MASK, (int2)(ind.x, ind.y)).w;
    #endif
    #endif
    #else
        maskBP[ind];
    #endif
}
#endif

#if defined(MASKFP)
DEVICE int readMaskFP(MASKFPTYPE maskFP, typeT ind) {
    return
#ifdef USEIMAGES
#ifdef CUDA
#ifdef MASKFP3D
        tex3D<unsigned char>(maskFP, ind.x, ind.y, ind.z);
#else
        tex2D<unsigned char>(maskFP, ind.x, ind.y);
#endif
#else
#ifdef MASKFP3D
        read_imageui(maskFP, sampler_MASK, (int4)(ind.x, ind.y, ind.z, 0)).w;
#else
        read_imageui(maskFP, sampler_MASK, (int2)(ind.x, ind.y)).w;
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
// #pragma OPENCL EXTENSION cl_ext_float_atomics : enable
// #define atomicAdd(a,b) atomic_fetch_add((volatile atomic_float *)(a),(b)) 
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

#ifdef TOF //////////////// TOF ////////////////
#define _2PI 0.3989423f

DEVICE float normPDF(const float x, const float mu, const float sigma) {

	const float a = (x - mu) / sigma;

	return _2PI / sigma * EXP(-0.5f * a * a);
}

DEVICE void TOFDis(const float3 diff, const float tc, const float LL, float* D, float* DD) {
	*D = length(diff * tc) - LL / 2.f;
	*DD = *D;
}

DEVICE float TOFWeight(const float element, const float sigma_x, const float D, const float DD, const float TOFCenter, float dX) {
	float output = normPDF(D, TOFCenter, sigma_x);
	dX *= sign(DD);
#pragma unroll
	for (int tr = 1; tr < CINT(TRAPZ_BINS) - 1; tr++)
#ifdef USEMAD
		output += (normPDF(FMAD(-dX, CFLOAT(tr), D), TOFCenter, sigma_x) * 2.f);
	output += normPDF(FMAD(-element, sign(DD), D), TOFCenter, sigma_x);
#else
		output += (normPDF(D - dX * CFLOAT(tr), TOFCenter, sigma_x) * 2.f);
	output += normPDF(D - element * sign(DD), TOFCenter, sigma_x);
#endif
	return output;
}


DEVICE float TOFLoop(const float DD, const float element, CONSTANT float* TOFCenter, const float sigma_x, float* D, const float epps, float* TOFWeights) {
	float TOFSum = 0.f;
	const float dX = element / (TRAPZ_BINS - 1.f);
#ifndef __CUDACC__ 
#pragma unroll NBINS
#endif
	for (int to = 0; to < NBINS; to++) {
		TOFWeights[to] = TOFWeight(element, sigma_x, *D, DD, TOFCenter[to], dX) * dX;
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
#ifdef CUDA
#ifdef USEIMAGES
    // if (local_ind.x <= 1.f && local_ind.y <= 1.f && local_ind.z <= 1.f && local_ind.x >= 0.f && local_ind.y >= 0.f && local_ind.z >= 0.f)
		*ax = (local_ele * tex3D<float>(d_OSEM, local_ind.x, local_ind.y, local_ind.z));
#else
	*ax = (local_ele * d_OSEM[local_ind]);
#endif
#else
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
#endif
}

// Computes the forward projection
// Includes TOF-specific weighting
DEVICE void denominator(PTR_THR float *ax, const typeT localInd, float local_ele, IMTYPE d_OSEM
#ifdef TOF
	, const float element, const float TOFSum, const float DD, const float sigma_x, float* D, float* TOFWeights
#ifdef LISTMODE
	, const int TOFIndex
#endif
#endif
#ifdef N_RAYS
	, const int lor
#endif
) {
	float apu = 0.f;
	forwardProject(local_ele, &apu, localInd, d_OSEM);
#ifdef TOF
	const float dX = apu / TOFSum;
#if defined(LISTMODE) && !defined(SENS)
	int to = TOFIndex;
#else
#ifndef __CUDACC__ 
#pragma unroll NBINS
#endif
	for (int to = 0; to < NBINS; to++) {
#endif
		const float joku = TOFWeights[to] * dX;
#ifdef N_RAYS
		ax[to + NBINS * lor] += joku;
#else
		ax[to] += joku;
#endif
#if !defined(LISTMODE) || defined(SENS)
	}
#endif
#else
#ifdef N_RAYS
	ax[lor] += apu;
#else
	ax[0] += apu;
#endif
#endif
}
#endif

#if defined(BP) && !defined(PROJ5) && (defined(ATOMIC) || defined(ATOMIC32) || defined(ATOMICF))
// Compute the backprojection
DEVICE void rhs(const float local_ele, PTR_THR const float *ax, const LONG local_ind, CLGLOBAL CAST* d_rhs_OSEM, const uchar no_norm, CLGLOBAL CAST* d_Summ
#ifdef TOF
	, const FLOAT element, const FLOAT sigma_x, float* D, const FLOAT DD, const FLOAT TOFSum, float* TOFWeights
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
#ifndef __CUDACC__ 
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
	atomic_add(&d_rhs_OSEM[local_ind], CINT(yaxTOF * TH));
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
		atomic_add(&d_Summ[local_ind], CINT(val * TH));
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
	const CLGLOBAL float* d_xy,
#else
	CONSTANT float* d_xy, 
#endif
	CONSTANT float* d_z, const CLGLOBAL ushort* trIndex, const CLGLOBAL ushort* axIndex, float3* s, float3* d, const size_t idx
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
	CONSTANT float* d_uv, float3* s, float3* d, const int3 i, const uint d_size_x, const uint d_sizey, const float2 d_dPitch
#ifdef PROJ5
	, float3* dR, float3* dL, float3* dU, float3* dD
#endif
) {
	int id = i.z * 6;
	*s = CMFLOAT3(d_xyz[id], d_xyz[id + 1], d_xyz[id + 2]);
	*d = CMFLOAT3(d_xyz[id + 3], d_xyz[id + 4], d_xyz[id + 5]);
	const float2 indeksi = MFLOAT2(CFLOAT(i.x) - CFLOAT(d_size_x) / 2.f + .5f, CFLOAT(i.y) - CFLOAT(d_sizey) / 2.f + .5f);
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
}

#elif defined(SPECT)
DEVICE void getDetectorCoordinatesSPECT(
#if defined(USEGLOBAL)
	const CLGLOBAL float *d_xyz,
#else
	CONSTANT float *d_xyz, 
#endif
    CONSTANT float *d_uv, PTR_THR FLOAT3 *s, PTR_THR FLOAT3 *d, const int3 i, const uint d_size_x, const uint d_sizey, const FLOAT2 d_dPitch, const CLGLOBAL float* d_rayShiftsDetector, const CLGLOBAL float* d_rayShiftsSource, int lorXY, size_t idx) {
	uint id = i.z * 6;
	*s = CMFLOAT3((FLOAT)d_xyz[id], (FLOAT)d_xyz[id + 1], (FLOAT)d_xyz[id + 2]); // TODO remove cast
	*d = CMFLOAT3((FLOAT)d_xyz[id + 3], (FLOAT)d_xyz[id + 4], (FLOAT)d_xyz[id + 5]); // TODO remove cast
	const FLOAT2 shift_det_elem = MFLOAT2(
        d_dPitch.x * (CFLOAT(i.x) + (FLOAT_ONE - CFLOAT(d_size_x)) * FLOAT_HALF),
        d_dPitch.y * (CFLOAT(i.y) + (FLOAT_ONE - CFLOAT(d_sizey)) * FLOAT_HALF)
    ); // Amount of shift from sinogram center to current detector element
	
    id = i.z * NA; // Index of d_uv (detector panel normal vector)
    uint idShift = 2*lorXY + (2*N_RAYS) * idx; // Index of rayShiftsDetector

	const FLOAT apuX = d_uv[id]; // X component of detector panel normal vector
	const FLOAT apuY = d_uv[id + 1]; // Y component of detector panel normal vector
    
	(*d).x += apuX * (shift_det_elem.x + d_rayShiftsDetector[idShift]); // Shift to current element + shift to rayShiftsDetector
	(*d).y += apuY * (shift_det_elem.x + d_rayShiftsDetector[idShift]);
	(*d).z += shift_det_elem.y + d_rayShiftsDetector[idShift+1];
	(*s).x += apuX * (shift_det_elem.x + d_rayShiftsSource[idShift]);
	(*s).y += apuY * (shift_det_elem.x + d_rayShiftsSource[idShift]);
	(*s).z += shift_det_elem.y + d_rayShiftsSource[idShift+1];

    *d += 100.f * (*d - *s); // Extend rays across FOV (100 too large for fp16 type)
}
#else
#if defined(RAW) || defined(SENS)
// Get the detector coordinates for the current (raw) measurement
DEVICE void getDetectorCoordinatesRaw(
#if defined(USEGLOBAL)
	const CLGLOBAL float* d_xy,
#else
	CONSTANT float* d_xy, 
#endif
	CONSTANT float *d_z, const int3 i, float3* s, float3* d, const int2 indz
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
#if defined(CT) || !defined(USEGLOBAL)
	CONSTANT float *d_xy, 
#else
	const CLGLOBAL float *d_xy, 
#endif
	CONSTANT float *d_z
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
	const CLGLOBAL float* d_xy,
#else
	CONSTANT float* d_xy, 
#endif
	CONSTANT float* d_z
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
#ifdef CUDA
#ifdef USEIMAGES
		tex3D<float>(d_atten, ind.x, ind.y, ind.z);
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
	int tempi = 0;
	float start = d_b - d + d_d;
	for (uint ii = 0u; ii < d_N; ii++) {
		if (start > 0.f) {
			tempi = CINT(ii);
			break;
		}
		start += d_d;
	}
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
			atnind += CLONG_rtz(iii * d_NN);
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
DEVICE float compute_element(PTR_THR float *t0, PTR_THR float *tc, const float L, const float tu, const int u, PTR_THR int *temp_ijk) {
	float local_ele = voxelValue(*t0, *tc, L);
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
