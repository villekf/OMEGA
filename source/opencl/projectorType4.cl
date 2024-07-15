
//CONSTANT sampler_t sampler_MASK = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_NONE;
//#ifdef SPECTMASK
//CONSTANT sampler_t sampler_SPECT = CLK_NORMALIZED_COORDS_TRUE | CLK_FILTER_NEAREST | CLK_ADDRESS_CLAMP;
//#endif

#ifndef NVOXELS
#define NVOXELS 1
#endif
//#define NVOXELS 8
//#define PITCH

#if defined(FP) || (defined(BP) && !defined(CT))// START FP
#define NSTEPS 10000
// #undef PITCH
// #undef NA
// #define NA 6

// Forward projection
KERNEL
void projectorType4Forward(const uint d_size_x, const uint d_sizey, 
#ifdef PYTHON
	const float d_dPitchX, const float d_dPitchY, 
#else
	const float2 d_dPitch, 
#endif
    const float dL, const float global_factor, 
#ifdef MASKFP
    IMAGE2D maskFP,
#endif
#if defined(MASKBP) && defined(BP) && !defined(CT)
    IMAGE2D maskBP,
#endif
	///////////////////////// TOF BINS /////////////////////////
#ifdef TOF
	CONSTANT float* TOFCenter, const float sigma_x, 
#endif
	///////////////////////// END TOF BINS /////////////////////////
    ////////////////////////////////////////////////////////////////////////
#if !defined(CT) && defined(ATN) && !defined(ATNM)
	IMAGE3D d_atten,
#elif !defined(CT) && !defined(ATN) && defined(ATNM)
	const CLGLOBAL float* CLRESTRICT d_atten,
#endif
    ////////////////////////////////////////////////////////////////////////
#ifdef PYTHON
	const uint d_Nx, const uint d_Ny, const uint d_Nz, const float bx, const float by, const float bz, 
    const float d_bmaxx, const float d_bmaxy, const float d_bmaxz, const float d_scalex, const float d_scaley, const float d_scalez,
#else
	const uint3 d_N, const float3 b, const float3 bmax, const float3 d_scale, 
#endif
#ifdef FP
    IMAGE3D d_OSEM, CLGLOBAL float* CLRESTRICT d_output,
#else
    const CLGLOBAL float* CLRESTRICT d_OSEM, CLGLOBAL CAST* CLRESTRICT d_output,
#endif
#if defined(CT) || defined(RAW) || defined(SUBSETS) || defined(SENS)
	CONSTANT float* d_xyz,
#else
	const CLGLOBAL float* CLRESTRICT d_xyz,
#endif
#if defined(LISTMODE) && !defined(SENS)
    const CLGLOBAL float* CLRESTRICT d_uv,
#else
    CONSTANT float* d_uv,
#endif
#ifdef SENS
	const int rings, const uint d_det_per_ring,
#endif
    const LONG d_nProjections,
#if defined(SUBSETS) && !defined(LISTMODE)
    const CLGLOBAL uint* CLRESTRICT d_xyindex, const CLGLOBAL ushort* CLRESTRICT d_zindex,
#endif
#ifdef RAW
    const CLGLOBAL ushort* CLRESTRICT d_L, const uint d_det_per_ring,
#endif
	///////////////////////// PET NORMALIZATION DATA /////////////////////////
#ifdef NORM
    const CLGLOBAL float* CLRESTRICT d_norm,
#endif
	///////////////////////// END PET NORMALIZATION DATA /////////////////////////
	///////////////////////// EXTRA CORRECTION DATA /////////////////////////
#ifdef SCATTER
    const CLGLOBAL float* CLRESTRICT d_scat,
#endif
	///////////////////////// END EXTRA CORRECTION DATA /////////////////////////
#if defined(BP) && !defined(CT)
    CLGLOBAL CAST* CLRESTRICT d_Summ,
#endif
	const uchar no_norm, const ULONG m_size, const uint currentSubset, const int aa)
{
    //
    int3 i = MINT3(GID0, GID1, GID2);
#if STYPE == 1 || STYPE == 2 || STYPE == 4 || STYPE == 5
    getIndex(&i, d_size_x, d_sizey, currentSubset);
#endif

#if (defined(CT) || defined(SPECT) || defined(PET)) && !defined(LISTMODE)
    size_t idx = GID0 + GID1 * d_size_x + GID2 * d_sizey * d_size_x;
    if (i.x >= d_size_x || i.y >= d_sizey || i.z >= d_nProjections)
#elif defined(SENS)
	int2 indz;
	indz.x = i.z / rings;
	indz.y = i.z % rings;
	size_t idx = i.x + i.y * d_det_per_ring + i.z * d_det_per_ring * d_det_per_ring;
	if (i.x >= d_det_per_ring || i.y >= d_det_per_ring || i.z >= rings * rings)
#elif defined(LISTMODE)
	size_t idx = GID0;
	if (idx >= m_size)
#else
    size_t idx = GID0 + GID1 * GSIZE0 + GID2 * GSIZE1 * GSIZE0;
    if (i.x >= d_size_x || idx >= m_size)
#endif
        return; 
#ifdef MASKFP
#ifdef CUDA
    const int maskVal = tex2D<unsigned char>(maskFP, i.x, i.y);
#else
    const int maskVal = read_imageui(maskFP, sampler_MASK, (int2)(i.x, i.y)).w;
#endif
    if (maskVal == 0)
        return;
#endif
#ifdef PYTHON
	const float2 d_dPitch = make_float2(d_dPitchX, d_dPitchY);
	const uint3 d_N = make_uint3(d_Nx, d_Ny, d_Nz);
	const float3 d_scale = make_float3(d_scalex, d_scaley, d_scalez);
	const float3 b = make_float3(bx, by, bz);
	const float3 bmax = make_float3(d_bmaxx, d_bmaxy, d_bmaxz);
#endif
#if defined(N_RAYS) && defined(FP)
	float ax[NBINS * N_RAYS];
#else
	float ax[NBINS];
#endif
#if defined(BP) && !defined(CT)
#ifndef __CUDACC__ 
#pragma unroll NBINS
#endif
	for (int to = 0; to < NBINS; to++)
#ifdef SENS
		ax[to] = 1.f;
#else
		ax[to] = d_OSEM[idx + to * m_size];
#endif
	// const float input = d_OSEM[idx];
#else
#if defined(N_RAYS) && defined(FP)
#ifndef __CUDACC__ 
#pragma unroll NBINS * N_RAYS
#endif
	for (int to = 0; to < NBINS * N_RAYS; to++)
		ax[to] = 0.f;
#else
#ifndef __CUDACC__ 
#pragma unroll NBINS
#endif
	for (int to = 0; to < NBINS; to++)
		ax[to] = 0.f;
#endif
#endif

#ifndef CT

#ifdef NORM // Normalization included
	float local_norm = 0.f;
	local_norm = d_norm[idx];
#endif
#ifdef SCATTER // Scatter data included
	float local_scat = 0.f;
	local_scat = d_scat[idx];
#endif
#endif

#ifdef N_RAYS //////////////// MULTIRAY ////////////////
	int lor = -1;
	// bool pass = true;
	// Load the next detector index
	// raw list-mode data
//#pragma unroll N_RAYS3D
	for (int lorZ = 0u; lorZ < N_RAYS3D; lorZ++) {
//#pragma unroll N_RAYS2D
		for (int lorXY = 0u; lorXY < N_RAYS2D; lorXY++) {
			lor++;
#endif  //////////////// END MULTIRAY ////////////////
	float3 s, d;
#if (defined(CT) || defined(SPECT)) && !defined(LISTMODE) && !defined(PET) // CT data
	getDetectorCoordinatesCT(d_xyz, d_uv, &s, &d, i, d_size_x, d_sizey, d_dPitch);
#elif defined(LISTMODE) && !defined(SENS) // Listmode data
	getDetectorCoordinatesListmode(d_xyz, &s, &d, idx
#if defined(N_RAYS)
		, lorXY, lorZ, d_dPitch
#endif
	);
#elif defined(RAW) || defined(SENS) // raw data
	getDetectorCoordinatesRaw(d_xyz, d_uv, i, &s, &d, indz
#if defined(N_RAYS)
		, lorXY, lorZ, d_dPitch
#endif
	);
#elif !defined(SUBSETS) // Precomputation phase, subset types 1, 2, 4, 5, 8, 9, 10, 11
	getDetectorCoordinatesFullSinogram(d_size_x, i, &s, &d, d_xyz, d_uv
#if defined(N_RAYS)
		, lorXY, lorZ, d_dPitch
#endif
#if defined(NLAYERS)
		, d_sizey
#endif
	);
#else // Subset types 3, 6, 7
	getDetectorCoordinates(d_xyindex, d_zindex, idx, &s, &d, d_xyz, d_uv
#if defined(N_RAYS)
		, lorXY, lorZ, d_dPitch
#endif
#if defined(NLAYERS)
		, d_sizey, d_size_x
#endif
	);
#endif
    float3 v = d - s;
    const float3 bmin = b;
    //const float3 bmax = CFLOAT(d_N) * d_d + b;
    const float3 tBack = (bmin - s) / v;
    const float3 tFront = (bmax - s) / v;
    //const float3 tBack = DIVIDE(bmin - s, v);
    //const float3 tFront = DIVIDE(bmax - s, v);

    const float3 tMin = fmin(tFront, tBack);
    const float3 tMax = fmax(tFront, tBack);

    const float tStart = fmax(fmax(tMin.x, tMin.y), tMin.z);
    const float tEnd = fmin(fmin(tMax.x, tMax.y), tMax.z);
#ifdef TOF //////////////// TOF ////////////////
    float TOFSum = 0.f;
#endif //////////////// END TOF ////////////////
    if (tStart >= tEnd)
#ifdef N_RAYS //////////////// MULTIRAY ////////////////
		continue;
#else
		return;
#endif  //////////////// END MULTIRAY ////////////////
	float L = length(v);
#ifndef CT
#ifdef ATN // Attenuation included
	float jelppi = 0.f;
#endif
#ifdef TOF //////////////// TOF ////////////////
	float D = 0.f;
	float DD = 0.f;
	TOFDis(v, tStart, L, &D, &DD);
#endif //////////////// END TOF ////////////////
#endif
    const float tStep = DIVIDE(dL, L);

    s = (s - bmin) * d_scale;
    v *= d_scale;
    float temp = 0.f;

    float t = tStart;
#if defined(ATN) && defined(BP)
    for (uint ii = 0; ii < NSTEPS; ii++) {
#ifndef USEMAD
        const float3 p = t * v + s;
#else
        const float3 p = FMAD(t, v, s);
#endif
		compute_attenuation(dL, p, d_atten, &jelppi, aa, d_N);
        t += tStep;
        if (t >= tEnd)
            break;
    }
    t = tStart;
#endif
#if !defined(CT) //////////////// PET ////////////////
#ifdef N_RAYS //////////////// MULTIRAY ////////////////
		temp = 1.f / (L * CFLOAT(N_RAYS));
#else
		temp = 1.f / L;
#endif //////////////// END MULTIRAY ////////////////
#if defined(ATN) && defined(BP)
		temp *= EXP(jelppi);
#endif
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
#endif //////////////// END PET ////////////////
    for (uint ii = 0; ii < NSTEPS; ii++) {
#ifndef USEMAD
        const float3 p = v * t + s;
#else
        const float3 p = FMAD3(t, v, s);
#endif

    // if (idx == 8000) {
        // printf("v.x = %f\n", v.x);
        // printf("v.y = %f\n", v.y);
        // printf("v.z = %f\n", v.z);
        // printf("s.x = %f\n", s.x);
        // printf("s.y = %f\n", s.y);
        // printf("s.z = %f\n", s.z);
        // printf("t = %f\n", t);
    // }
#if defined(ATN) && defined(FP)
		compute_attenuation(dL, p, d_atten, &jelppi, aa, d_N);
#endif
#ifdef TOF //////////////// TOF ////////////////
			TOFSum = TOFLoop(DD, dL, TOFCenter, sigma_x, &D, 1e-6f);
#endif //////////////// END TOF ////////////////
#if defined(FP) //////////////// FORWARD PROJECTION ////////////////
			denominator(ax, p, dL, d_OSEM
#ifdef TOF //////////////// TOF ////////////////
			, dL, TOFSum, DD, TOFCenter, sigma_x, &D
#endif //////////////// END TOF ////////////////
#ifdef N_RAYS
			, lor
#endif
			);
    // if (idx == 8000) {
        // printf("ax[0] = %f\n", ax[0]);
        // printf("dL = %f\n", dL);
        // printf("p.x = %f\n", p.x);
        // printf("p.y = %f\n", p.y);
        // printf("p.z = %f\n", p.z);
    // }
#endif  //////////////// END FORWARD PROJECTION ////////////////
#if defined(BP) && !defined(CT) //////////////// BACKWARD PROJECTION ////////////////
            const uint local_ind = CUINT_rtz(p.x * CFLOAT(d_N.x)) + CUINT_rtz(p.y * CFLOAT(d_N.y)) * d_N.x + CUINT_rtz(p.z * CFLOAT(d_N.z)) * d_N.x * d_N.y;
            if (local_ind <= 0 || local_ind >= d_N.x * d_N.y * d_N.z) {
                t += tStep;
                if (t > tEnd)
                    break;
                else
                    continue;
            }
#if defined(MASKBP)
            int maskVal = 1;
            if (aa == 0) {
#ifdef CUDA
		        maskVal = tex2D<unsigned char>(maskBP, p.x, p.y);
#else
		        maskVal = read_imageui(maskBP, sampler_MASK, (int2)(p.x, p.y)).w;
#endif
            }
            if (maskVal > 0)
#endif
			rhs(dL * temp, ax, local_ind, d_output, no_norm, d_Summ
#ifdef TOF
			, dL, sigma_x, &D, DD, TOFCenter, TOFSum
#endif
			);
#endif //////////////// END BACKWARD PROJECTION ////////////////
        t += tStep;
        if (t > tEnd)
            break;
    }
    // if (idx == 8000) {
    //     printf("ax[0] = %f\n", ax[0]);
    //     printf("d_output[idx] = %f\n", d_output[idx]);
    // }
#if defined(ATN) && defined(FP)
		temp *= EXP(jelppi);
#endif
#if defined(FP) && !defined(N_RAYS) //////////////// FORWARD PROJECTION ////////////////
#ifndef __CUDACC__ 
#pragma unroll NBINS
#endif
		for (size_t to = 0; to < NBINS; to++) {
			forwardProjectAF(d_output, ax, idx, temp, to);
#ifdef TOF
			idx += m_size;
#endif
		}
    // if (idx == 8000) {
    //     printf("d_output[idx] = %f\n", d_output[idx]);
    // }
#elif defined(FP) && defined(N_RAYS)
#ifndef __CUDACC__ 
#pragma unroll NBINS
#endif
	for (int to = 0; to < NBINS; to++)
		ax[to + NBINS * lor] *= temp;
#endif //////////////// END FORWARD PROJECTION ////////////////
#ifdef N_RAYS //////////////// MULTIRAY ////////////////
		}
	}
#if defined(FP) //////////////// FORWARD PROJECTION ////////////////
#ifndef __CUDACC__ 
#pragma unroll NBINS
#endif
    for (size_t to = 0; to < NBINS; to++) {
        float apu = 0.f;
#ifndef __CUDACC__ 
#pragma unroll N_RAYS
#endif
        for (size_t kk = 0; kk < N_RAYS; kk++) {
            apu += ax[to + NBINS * kk];
        }
        ax[to] = apu;
    }
#ifndef __CUDACC__ 
#pragma unroll NBINS
#endif
    for (size_t to = 0; to < NBINS; to++) {
        forwardProjectAF(d_output, ax, idx, 1.f, to);
#ifdef TOF
        idx += m_size;
#endif
    }
#endif //////////////// END FORWARD PROJECTION ////////////////
#endif //////////////// END MULTIRAY ////////////////
    // if (idx == 8000) {
    //     printf("d_output[idx] = %f\n", d_output[idx]);
    // }
}
#endif // END FP

#if defined(BP) && defined(CT)// START BP
//#define MAX(a,b) (a>b?a:b)
//#define MIN(a,b) (a<b?a:b)
// #undef PITCH
// #undef NA
// #define NA 6
KERNEL2
void projectorType4Backward(const uint d_size_x, const uint d_sizey, 
#ifdef PYTHON
	const float d_dPitchX, const float d_dPitchY, 
#else
	const float2 d_dPitch, 
#endif
#ifdef USEIMAGES
#ifdef MASKBP
    IMAGE2D maskBP,
#endif
#else
#ifdef MASKBP
    const CLGLOBAL uchar* CLRESTRICT maskBP,
#endif
#endif
#ifdef OFFSET
    // const float T, //const float R, 
    CONSTANT float* T, //const float R, 
//     IMAGE2D maskOffset,
#endif
#ifdef PYTHON
	const uint d_Nx, const uint d_Ny, const uint d_Nz, const float bx, const float by, const float bz, 
    const float d_dx, const float d_dy, const float d_dz,
#else
	const uint3 d_N, const float3 b, const float3 d_d, 
#endif
    const float kerroin, 
#ifdef USEIMAGES
    IMAGE3D d_forw, 
#else
    const CLGLOBAL float* CLRESTRICT d_forw, 
#endif
#ifdef FDK
    CONSTANT float* angle, const float DSC, 
#endif
    CLGLOBAL float* CLRESTRICT d_OSEM, CONSTANT float* d_xyz, CONSTANT float* d_uv, CLGLOBAL float* CLRESTRICT d_Summ, 
    const uchar no_norm, const LONG d_nProjections, const int ii) {

    const int3 i = MINT3(GID0, GID1, GID2 * NVOXELS);

#ifdef PYTHON
	const uint3 d_N = make_uint3(d_Nx, d_Ny, d_Nz);
#endif
    if (i.x >= d_N.x || i.y >= d_N.y || i.z >= d_N.z)
        return;
    size_t idx = GID0 + GID1 * d_N.x + GID2 * NVOXELS * d_N.y * d_N.x;
#ifdef MASKBP
    if (ii == 0) {
#ifdef USEIMAGES
#ifdef CUDA
        const int maskVal = tex2D<unsigned char>(maskBP, i.x, i.y);
#else
        const int maskVal = read_imageui(maskBP, sampler_MASK, (int2)(i.x, i.y)).w;
#endif
#else
        const int maskVal = maskBP[i.x + i.y * d_N.x];
#endif
        if (maskVal == 0)
            return;
    }
#endif
#ifdef PYTHON
	const float2 d_dPitch = make_float2(d_dPitchX, d_dPitchY);
	const float3 d_d = make_float3(d_dx, d_dy, d_dz);
	const float3 b = make_float3(bx, by, bz);
#endif
    float temp[NVOXELS];
    float wSum[NVOXELS];
// #ifdef OFFSET
//     __private float tempOff[NVOXELS];
//     __private float wSumOff[NVOXELS];
// #endif
    for (int zz = 0; zz < NVOXELS; zz++) {
        temp[zz] = 0.f;
        if (no_norm == 0u)
            wSum[zz] = 0.f;
    }
// #ifdef OFFSET
//     __private float TT[600] = {11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.398f,11.12f,11.12f,11.398f,11.12f,11.398f,11.398f,11.12f,11.398f,11.398f,11.12f,11.398f,11.12f,11.398f,11.398f,11.12f,11.398f,11.398f,11.12f,11.398f,11.398f,11.398f,11.398f,11.12f,11.398f,11.398f,11.12f,11.398f,11.398f,11.12f,11.398f,11.12f,11.398f,11.398f,11.12f,11.12f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.398f,11.12f,11.398f,11.398f,11.12f,11.398f,11.398f,11.398f,11.398f,11.12f,11.398f,11.398f,11.12f,11.398f,11.398f,11.12f,11.398f,11.398f,11.398f,11.398f,11.12f,11.12f,11.12f,11.12f,11.398f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.398f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,10.842f,11.12f,11.12f,10.842f,11.12f,10.842f,10.842f,11.12f,10.842f,10.842f,11.12f,10.842f,11.12f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,11.12f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,10.842f,11.12f,11.12f,10.842f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,10.842f,11.12f,10.842f,11.12f,11.12f,10.842f,11.12f,11.12f,11.12f,11.12f,10.842f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f,11.12f};
// #endif
    //const float nZ = 1.f / CFLOAT(d_nProjections);
    //float pz = nZ / 2.f;
    float3 dV = CFLOAT3(i) * d_d + d_d / 2.f + b;
    const float2 koko = MFLOAT2(CFLOAT(d_size_x) * d_dPitch.x, CFLOAT(d_sizey) * d_dPitch.y );
    const float2 indeksi = MFLOAT2(CFLOAT(d_size_x) / 2.f, CFLOAT(d_sizey) / 2.f );
    //float apuZ = d_d.z / 2.f + b.z;
//#pragma unroll NPROJECTIONS
    for (int kk = 0; kk < d_nProjections; kk++) {
        float3 d1, d2, d3;
        float3 s;
        //const float kulmaXY = d_angles[kk * NA];
        s = CMFLOAT3(d_xyz[kk * 6], d_xyz[kk * 6 + 1], d_xyz[kk * 6 + 2]);
        d1 = CMFLOAT3(d_xyz[kk * 6 + 3], d_xyz[kk * 6 + 4], d_xyz[kk * 6 + 5]);
// #ifdef OFFSET
//         const float R = distance(s, d1);
// #endif
        //s = vload3(kk * 2, d_xyz);

        //d1 = vload3(kk * 2, &d_xyz[3]);
#if defined(PITCH)
        const float3 apuX = CMFLOAT3(d_uv[kk * NA], d_uv[kk * NA + 1], d_uv[kk * NA + 2]) * indeksi.x;
        const float3 apuY = CMFLOAT3(d_uv[kk * NA + 3], d_uv[kk * NA + 4], d_uv[kk * NA + 5]) * indeksi.y;
        //d2.x = d1.x + indeksi.x * d_uv[kk * NA] + indeksi.y * d_uv[kk * NA + 1];
        //d3.x = d1.x - indeksi.x * d_uv[kk * NA] - indeksi.y * d_uv[kk * NA + 1];
        //d2.y = d1.y + indeksi.x * d_uv[kk * NA + 2] + indeksi.y * d_uv[kk * NA + 3];
        //d3.y = d1.y - indeksi.x * d_uv[kk * NA + 2] - indeksi.y * d_uv[kk * NA + 3];
        //d2.z = d1.z + indeksi.x * d_uv[kk * NA + 4] + indeksi.y * d_uv[kk * NA + 5];
        //d3.z = d1.z - indeksi.x * d_uv[kk * NA + 4] - indeksi.y * d_uv[kk * NA + 5];
#else
        const float3 apuX = MFLOAT3(indeksi.x * d_uv[kk * NA], indeksi.x * d_uv[kk * NA + 1], 0.f);
        const float3 apuY = MFLOAT3(0.f, 0.f, indeksi.y * d_dPitch.y);
#endif
        d2 = apuX - apuY;
        d3 = d1 - apuX - apuY;
        //const float3 normX = normalize(apuX);
        //const float3 normY = normalize(apuY);
        const float3 normX = normalize(apuX) / koko.x;
        const float3 normY = normalize(apuY) / koko.y;
        const float3 cP = cross(d2, d3 - d1);
        const float upperPart = dot(cP, s - d1);
        float3 v = dV - s;
        float lowerPart = -dot(v, cP);
        //const float3 dMin = d2 - d3;
#ifndef USEMAD
        const float vApu = v.x * v.x + v.y * v.y;
#else
        const float vApu = FMAD(v.x, v.x, v.y * v.y);
#endif
        const float dApu = d_d.z * cP.z;
        const float pz = (CFLOAT(kk) + 0.5f) / CFLOAT(d_nProjections);
#ifndef __CUDACC__ 
#pragma unroll NVOXELS
#endif
        for (int zz = 0; zz < NVOXELS; zz++) {
            //d.z = dV.z;
            const uint ind = i.z + zz;
            if (ind >= d_N.z)
                break;
            const float t = DIVIDE(upperPart, lowerPart);
#ifndef USEMAD
            float3 p = s + v * t;
            const float l1 = vApu + v.z * v.z;
#else
            float3 p = FMAD3(t, v, s);
            const float l1 = FMAD(v.z, v.z, vApu);
#endif
#ifdef FDK
            // float weight = (DSC + dV.x * SINF(angle[kk]) - dV.y * COSF(angle[kk]));
            float weight = (DSC + dV.x * COSF(angle[kk]) - dV.y * SINF(angle[kk]));
            weight = (DSC * DSC) / (weight * weight) * (M_PI_F / (CFLOAT(d_nProjections) * d_dPitch.x));
            // weight = (DSC * DSC) / (weight * weight) * (2.f * M_PI_F / (CFLOAT(d_nProjections) * d_dPitch.x));
            // const float L = distance(p, s);
            // const float weight = (L * L * L) / (l1)*kerroin / 6.f;
#else
            const float L = distance(p, s);
            //const float l1 = v.x * v.x + v.y * v.y + v.z * v.z;
            const float weight = (L * L * L) / (l1)*kerroin;
#endif
            p -= d3;
            //const float yVar = read_imagef(d_forw, samplerIm, CFLOAT4(dot(p, normX) / koko.x, dot(p, normY) / koko.y, pz, 0.f)).x;
            float px = dot(p, normX);
            float py = dot(p, normY);
            float yVar = 0.f;
#ifdef USEIMAGES
#ifdef CUDA
            if (px <= 1.f && py <= 1.f && pz <= 1.f && px >= 0.f && py >= 0.f && pz >= 0.f)
                yVar = tex3D<float>(d_forw, px, py, pz);
#else
            if (px <= 1.f && py <= 1.f && pz <= 1.f && px >= 0.f && py >= 0.f && pz >= 0.f)
                yVar = read_imagef(d_forw, samplerIm, CFLOAT4(px, py, pz, 0.f)).w;
#endif
#else
            if (px < 1.f && py < 1.f && pz < 1.f && px >= 0.f && py >= 0.f && pz >= 0.f) {
                const LONG indX = CLONG_rtz(px * CFLOAT(d_size_x));
                const LONG indY = CLONG_rtz(py * CFLOAT(d_sizey)) * CLONG_rtz(d_size_x);
                const LONG indZ = CLONG_rtz(pz * CFLOAT(d_nProjections)) * CLONG_rtz(d_sizey) * CLONG_rtz(d_size_x);
                yVar = d_forw[indX + indY + indZ];
            }
#endif
			// if (idx == 25824640) {
				// printf("yVar = %f\n", yVar);
				// printf("px = %f\n", px);
				// printf("py = %f\n", py);
				// printf("weight = %f\n", weight);
				// printf("p.x = %f\n", p.x);
				// printf("p.y = %f\n", p.y);
				// printf("p.z = %f\n", p.z);
			// }
#ifdef OFFSET
            // const int offsetMask = read_imageui(maskOffset, sampler_Offset, (float2)(px, py)).w;
            // const int offsetMask = 1;
            // const float TT = 0.278f * 50.f;
            // const float R = 1.083999e3f;
            // const float apx = px;
            float TT;
            float Tloc = T[kk];
            if (Tloc > koko.x / 2.f) {
                px = fabs(px - 1.f);
                TT = koko.x - Tloc;
            }
            else
                TT = Tloc;
            px *= koko.x;
            px -= TT;
            // px -= TT[kk];
            // p -= T;
#endif
                // if (ii == 3 && i.x == 15 && i.y == 80 && i.z == 112) {
                //     printf("yVar = %f\n", yVar);
                //     printf("weight = %f\n", weight);
                // }
            //const float yVar = read_imagef(d_forw, sampler, CFLOAT4(p.x, p.z, pz, 0.f)).x;
            //const float yVar = p.x;
            if (yVar != 0.f) {
//                 float LO = 63.f;
//                 py *= koko.y;
//                 py -= koko.y / 2.f;
//                 py += LO;
//                 LO = koko.y / 2.f - LO;
// #ifdef MASKBP
                // if (maskVal > 0)
                // temp[zz] += yVar * weight *.5f;
                // else
// #endif
//                 float W2 = 1.f;
//                 if (py <= LO || py >= -LO)
//                     W2 = cos(M_PI_4_F * (atan(py / 700.f) / atan(LO / 700.f) - 1.f));
                // if (ii == 3 && i.x == 15 && i.y == 80 && i.z == 13) {
                //     // printf("dot(p, normX) = %f\n", dot(p, normX));
                //     // printf("dot(p, normY) = %f\n", dot(p, normY));
                //     printf("yVar = %f\n", yVar);
                //     printf("weight = %f\n", weight);
                //     // printf("px = %f\n", px);
                //     // printf("W2 = %f\n", W2);
                //     // printf("LO = %f\n", LO);
                //     // printf("koko.x = %f\n", koko.x);
                //     // printf("koko.y = %f\n", koko.y);
                //     // printf("kk = %d\n", kk);
                //     // printf("zz = %d\n", zz);
                // }
                // temp[zz] += yVar * weight * W2 * W2;
#ifdef OFFSET
                // if (ii < 3 && px <= TT[kk] && px >= -TT[kk]) {
                if (px <= TT && px >= -TT) {
                // if (px <= TT + 18.648f - 273.5040f && px >= TT - 273.5040f) {
                // if (px <= 273.5040f - 200.f && px >= -TT + 273.5040f) {
                    // px += TT;
                    // const float w = .5f * (SINF((M_PI_F * atan(px / R)) / (2.f * atan(TT[kk] / R))) + 1.f);
                    // const float w = .5f * (SINF((M_PI_F * atan(px / R)) / (2.f * atan(TT / R))) + 1.f);
                    // const float w = (.5f * (SINF((M_PI_F * atan2(px, R)) / (2.f * atan2(TT, R))) + 1.f));
                    // float w = SINF(M_PI_2_F * ((px + TT) / (2.f * TT)));
                    // w *= w;
                    float w = .5f * (1.f + SINF(M_PI_F * px / (TT * 2.f)));
                    // w /= 2.f;
                    // float w = .5f;
                    // float w = 10.f;
                    // temp[zz] += w * weight;
                    temp[zz] += w * yVar * weight;
                    //temp += yVar;
                    if (no_norm == 0u)
                        wSum[zz] += w * weight;
                // if (ind == 262 && i.y < 500 && i.y > 400 && i.x < 500 && i.x > 400) {
                //     printf("px = %f\n", px);
                //     printf("apx = %f\n", apx);
                //     printf("koko.x = %f\n", koko.x);
                //     printf("T = %f\n", T);
                //     printf("R = %f\n", R);
                //     printf("w = %f\n", w);
                //     printf("i.x = %d\n", i.x);
                //     printf("i.y = %d\n", i.y);
                // }
                }
                // else if (ii < 3 && px < -TT[kk]) {
                else if (px < -TT) {
                // else if (px < TT - 273.5040f) {
                // else if (px < -TT + 273.5040f) {
                }
                else {
#endif
                    // temp[zz] += weight;
                    temp[zz] += yVar * weight;
                    if (no_norm == 0u)
                        wSum[zz] += weight;
                    // if (i.x == 700 && i.y == 700 && ind == 50)  {
                    //    printf("px = %f\n", px);
                    //    printf("py = %f\n", py);
                    //    printf("pz = %f\n", pz);
                    //    printf("weight = %f\n", weight);
                    //    printf("yVar = %f\n", yVar);
                    //    printf("temp[zz] = %f\n", temp[zz]);
                    // }
#ifdef OFFSET
                }
#endif
            }
            v.z += d_d.z;
            //lowerPart -= d_d.z * cP.z;
            lowerPart -= dApu;
        }
        //pz += nZ;
        //break;
//#endif
    }
// #ifdef BP
    for (int zz = 0; zz < NVOXELS; zz++) {
        const uint ind = i.z + zz;
        if (ind >= d_N.z)
            break;
                // if (ii == 3 && i.x == 15 && i.y == 80 && i.z == 112) {
                //     printf("temp[zz] = %f\n", temp[zz]);
                //     printf("zz = %d\n", zz);
                // }
// #ifdef OFFSET
//         if (temp[zz] != 0.f) {
// #endif
            d_OSEM[idx] += temp[zz];

            if (no_norm == 0u)
                d_Summ[idx] = wSum[zz];
// #ifdef OFFSET
//         }
//         else {
//             d_OSEM[idx] = tempOff[zz];
//             if (no_norm == 0u)
//                 d_Summ[idx] = wSumOff[zz];
//         }
// #endif
        idx += d_N.y * d_N.x;
    }
// #endif
}
#endif // END BP