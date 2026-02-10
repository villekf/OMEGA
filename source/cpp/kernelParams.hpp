// Structs to save buffer slots with Metal kernels and scalar inputs. This is a hybrid C++/MSL/OpenCL/CUDA file with preprocessor directives. TODO: use with OpenCL and CUDA

#if defined(__METAL_VERSION__) // MSL
#define CL_FLOAT2 float2
#define CL_FLOAT3 float3
#define CL_INT3 int3
#define CL_INT32 int
#define CL_INT64 long
#define CL_UINT3 uint3
#define CL_UINT8 unsigned char
#define CL_UINT32 uint
#define CL_UINT64 unsigned long
#elif defined(__OPENCL_VERSION__) // OpenCL
#define CL_FLOAT2 float2
#define CL_FLOAT3 float3
#define CL_INT3 int3
#define CL_INT32 int
#define CL_INT64 long
#define CL_UINT3 uint3
#define CL_UINT8 unsigned char
#define CL_UINT32 uint
#define CL_UINT64 unsigned long
#elif defined(__CUDACC__) // CUDA
#elif (defined(METAL) && !defined(__METAL_VERSION__)) // C++ on MacOS
#define CL_FLOAT2 simd::float2
#define CL_FLOAT3 simd::float3
#define CL_INT3 simd::int3
#define CL_INT32 int
#define CL_INT64 long
#define CL_UINT3 simd::uint3
#define CL_UINT8 unsigned char
#define CL_UINT32 uint
#define CL_UINT64 unsigned long
#elif (defined(OPENCL) && !defined(__OPENCL_VERSION__)) // C++ / OpenCL
#define CL_FLOAT2 cl_float2
#define CL_FLOAT3 cl_float3
#define CL_INT3 cl_int3
#define CL_INT32 cl_int
#define CL_INT64 cl_long
#define CL_UINT3 cl_uint3
#define CL_UINT8 cl_uchar
#define CL_UINT32 cl_uint
#define CL_UINT64 cl_ulong
#elif (defined(CUDA) && !defined(__CUDACC__)) // C++ / CUDA
#endif

struct ScalarKernelParams { // Kernel scalar values that do not change with time step or (sub)iteration. See initializeKernel in ProjectorClass.h
    // Static
    CL_UINT32 nRowsD;
    CL_UINT32 nColsD;
    CL_FLOAT2 dPitch;
    float dL;
    float global_factor;
    float epps;
    CL_UINT32 det_per_ring;
    float sigma_x;
    float coneOfResponseStdCoeffA;
    float coneOfResponseStdCoeffB;
    float coneOfResponseStdCoeffC;
    float tube_width;
    float cylRadiusProj3;
    float bmin;
    float bmax;
    float Vmax;
    CL_UINT32 rings;
    float helicalRadius;
    // Dynamic (change per subset or timestep)
    CL_UINT3 d_N;
    CL_FLOAT3 b;
    CL_FLOAT2 dSize5;
    float kerroin4;
    float DSC;
    CL_FLOAT3 d;
    CL_FLOAT3 d_Scale4;
    CL_FLOAT3 d_Scale5;
    CL_FLOAT3 d_bmax;
    float orthWidth;
    CL_INT64 nProjections;
    CL_UINT8 no_norm;
    CL_UINT64 m_size;
    CL_UINT32 currentSubset;
    CL_INT32 aa;
};

#undef CL_FLOAT2
#undef CL_FLOAT3
#undef CL_INT3
#undef CL_INT32
#undef CL_INT64
#undef CL_UINT3
#undef CL_UINT8
#undef CL_UINT32
#undef CL_UINT64