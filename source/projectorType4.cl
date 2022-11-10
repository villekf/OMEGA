
//__constant sampler_t sampler_MASK = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_NONE;
//#ifdef SPECTMASK
//__constant sampler_t sampler_SPECT = CLK_NORMALIZED_COORDS_TRUE | CLK_FILTER_NEAREST | CLK_ADDRESS_CLAMP;
//#endif

//#define NVOXELS 8
//#define PITCH

#if defined(FP)
#define NSTEPS 10000

// Forward projection
__kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, 1)))
void projectorType4Forward(const uint3 d_N, const float3 b, const uint d_size_x, const uint d_sizey, const float2 d_dPitch, const float3 bmax, 
    const float dL, const float3 d_scale, const float global_factor,
#ifdef MASKFP
    __read_only image2d_t maskFP
#endif
    ////////////////////////////////////////////////////////////////////////
#if !defined(CT) && defined(ATN)
    __read_only image3d_t d_atten,
#endif
//    ////////////////////////////////////////////////////////////////////////
    __read_only image3d_t d_OSEM, __global float* restrict d_forw,
#ifdef CT
    __constant float* d_xyz,
#else
    const __global float* restrict d_xyz,
#endif
#if defined(LISTMODE)
    const __global float* restrict d_uv,
#else
    __constant float* d_uv,
#endif
    const long d_nProjections
#if defined(SUBSETS) && !defined(LISTMODE)
    , const __global uint* restrict d_xyindex, const __global ushort* restrict d_zindex
#endif
#ifdef RAW
    , const __global ushort* restrict d_L, const uint d_det_per_ring
#endif
    ////////////////////////////////////////////////////////////////////////
#ifdef NORM
    , const __global float* restrict d_norm
#endif
    ////////////////////////////////////////////////////////////////////////
#ifdef SCATTER
    , const __global float* restrict d_scat
#endif
    ////////////////////////////////////////////////////////////////////////
)
{
    ///*
    const int3 i = { get_global_id(0), get_global_id(1), get_global_id(2) };

#if defined(CT) || defined(SPECT) || defined(PET)
    size_t idx = get_global_id(0) + get_global_id(1) * d_size_x + get_global_id(2) * d_sizey * d_size_x;
    if (i.x >= d_size_x || i.y >= d_sizey || i.z >= d_nProjections)
#else
    size_t idx = get_global_id(0) + get_global_id(1) * get_global_size(0) + get_global_id(2) * get_global_size(1) * get_global_size(0);
    if (i.x >= d_size_x || idx >= m_size)
#endif
        return; 
#ifdef MASKFP
    const int maskVal = read_imageui(maskFP, sampler_MASK, (int2)(i.x, i.y)).x;
    if (maskVal == 0)
        return;
#endif
    float3 s, d;
#if !defined(CT) // PET and SPECT ONLY
    float d_epsilon_mramla = 0.f;
    float d_epps = 0.f;
    float local_sino = 0.f;
    const uint d_NN = d_N.x * d_N.y * d_N.z;
#endif
#ifdef LISTMODE
    getDetectorCoordinatesListmode(d_xyz, &s, &d, idx);
#else
#if defined(CT) //|| (defined(SPECT) && !defined(SPECTMASK))
    getDetectorCoordinatesCT(d_xyz, d_uv, &s, &d, i, d_size_x, d_sizey, d_dPitch);
#elif defined(RAW) // raw data
    getDetectorCoordinatesRaw(d_xyz, d_uv, d_L, d_det_per_ring, idx, &s, &d); // Sinogram data
#elif defined(FIND_LORS) || !defined(SUBSETS) // Precomputation phase
    getDetectorCoordinatesFullSinogram(d_size_x, i, &s, &d, d_xyz, d_uv);
#else // Not the precomputation phase
    getDetectorCoordinates(d_xyindex, d_zindex, d_size_x, idx, &s, &d, d_xyz, d_uv);
#endif
#endif
#if !defined(CT)
#ifdef ATN // Attenuation included
    float jelppi = 0.f;
#endif
    float local_norm = 0.f;
    float local_scat = 0.f;
#ifdef NORM // Normalization included
    local_norm = d_norm[idx];
#endif
#ifdef SCATTER // Scatter data included
    local_scat = d_scat[idx];
#endif
#ifdef AF
    float ax[NROLLS];
#pragma unroll
    for (uint kk = 0; kk < NROLLS; kk++)
        ax[kk] = 0.f;
#else
#ifdef TOF
    float ax[NBINS];
#pragma unroll NBINS
    for (uint to = 0; to < NBINS; to++)
        ax[to] = 0.f;
#else
    float axOSEM = 0.f;
#endif
#endif
#ifdef TOF
    float D = 0.f, DD = 0.f;
#endif
#endif

    float3 v = d - s;
    const float3 bmin = b;
    //const float3 bmax = convert_float(d_N) * d_d + b;
    const float3 tBack = (bmin - s) / v;
    const float3 tFront = (bmax - s) / v;
    //const float3 tBack = native_divide(bmin - s, v);
    //const float3 tFront = native_divide(bmax - s, v);

    const float3 tMin = fmin(tFront, tBack);
    const float3 tMax = fmax(tFront, tBack);

    const float tStart = fmax(fmax(tMin.x, tMin.y), tMin.z);
    const float tEnd = fmin(fmin(tMax.x, tMax.y), tMax.z);
    //if (i.x == 102 && i.y == 80 && i.z == 10) {
    //    printf("tStart = %f\n", tStart);
    //    printf("tEnd = %f\n", tEnd);
    //    printf("tMin.x = %f\n", tMin.x);
    //    printf("tMin.y = %f\n", tMin.y);
    //    printf("tMin.z = %f\n", tMin.z);
    //    printf("tMax.x = %f\n", tMax.x);
    //    printf("tMax.y = %f\n", tMax.y);
    //    printf("tMax.z = %f\n", tMax.z);
    //    printf("s.x = %f\n", s.x);
    //    printf("s.y = %f\n", s.y);
    //    printf("s.z = %f\n", s.z);
    //    printf("d.x = %f\n", d.x);
    //    printf("d.y = %f\n", d.y);
    //    printf("d.z = %f\n", d.z);
    //}
    if (tStart >= tEnd)
        return;
    const float tStep = native_divide(dL, length(v));
#ifndef CT
    //const float L = length(v);
    //float jelppi = 0.f;
    //float dP = 0.f;
    float dP = tStep;
#endif

    s = (s - bmin) * d_scale;
    v *= d_scale;
    //float3 p = (s + tStart * v - bmin) * d_scale;
    //float3 p = s + tStart * v;
#ifdef FP
    float temp = 0.f;

    float t = tStart;
    for (uint ii = 0; ii < NSTEPS; ii++) {
        const float4 p = (float4)(mad(t, v, s), 0.f);
#if !defined(CT)
#ifdef ATN
        jelppi -= read_imagef(d_atten, samplerForw, p).x;
#endif
#if defined(TOF) && (defined(DEC) || defined(FP))
        const float TOFSum = TOFLoop(DD, dL, store_elements, TOFCenter, sigma_x, &D, ii * NBINS, d_epps);
#endif
#ifdef FP
#ifdef AF // Implementation 2

         denominator(dL, ax, p, d_NN, d_OSEM);

#else // Implementation 3
            denominator_multi(dL, &axOSEM, read_imagef(d_OSEM, samplerForw, p).x);

#endif
#endif
            // dP += dL;
#else
        temp += read_imagef(d_OSEM, samplerForw, p).x;
#endif
        t += tStep;
        if (t > tEnd)
            break;
        //p = (s + t * v - bmin) * d_scale;
        //p = s + t * v;
        //p = mad(t, v, s);
    }

#ifndef CT
#ifdef ATN
    dP *= native_exp(jelppi * dL);
#endif
#ifdef NORM
    dP *= local_norm;
#endif
#ifdef SCATTER
    dP *= d_scat[idx];
#endif
    dP *= global_factor;

#ifdef FP
    //if (local_sino != 0.f) {
//#ifdef TOF
//        nominatorTOF(MethodList, ax, d_Sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx, m_size, local_sino);
//#else
#ifdef AF
        nominator(ax, local_sino, d_epsilon_mramla, d_epps, temp, idx
#if defined(BP) && defined(FP)
            , MethodList
#if defined(RANDOMS)
            , d_sc_ra
#endif
#endif
        );
        forwardProjectAF(d_forw, ax, idx, d_NN);
#else
        nominator_multi(&axOSEM, local_sino, d_epps, temp, idx
#if defined(RANDOMS) && defined(BP) && defined(FP)
            , d_sc_ra
#endif
        );
        d_forw[idx] = axOSEM;
#endif
//#endif
    //}
#endif
#else
//#if defined(FP) && defined(BP) && !defined(CT)
//    d_forw[idx] = d_Sino[idx] / (temp * dL);
//#else
    d_forw[idx] = temp * dL;
//#endif
#endif
#endif
//*/
}

#endif

#if defined(BP)
//#define MAX(a,b) (a>b?a:b)
//#define MIN(a,b) (a<b?a:b)
//#undef PITCH
//#undef NA
//#define NA 6
__kernel __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, 1)))
void projectorType4Backward(const uint3 d_N, const float3 b, const uint d_size_x, const uint d_sizey, const float2 d_dPitch, const float3 d_d, const float kerroin,
#ifdef MASKBP
    __read_only image2d_t maskBP,
#endif
    __read_only image3d_t d_forw, __global float* restrict d_OSEM, __global float* restrict d_Summ, __constant float* d_xyz, __constant float* d_uv,
    const uchar no_norm, const long d_nProjections) {

    const uint3 i = { get_global_id(0), get_global_id(1), get_global_id(2) * NVOXELS };
    size_t idx = get_global_id(0) + get_global_id(1) * d_N.x + get_global_id(2) * NVOXELS * d_N.y * d_N.x;

    if (i.x >= d_N.x || i.y >= d_N.y || i.z >= d_N.z)
        return;
#ifdef MASKBP
    const int maskVal = read_imageui(maskBP, sampler_MASK, (int2)(i.x, i.y)).x;
    if (maskVal == 0)
        return;
#endif

    __private float temp[NVOXELS];
    __private float wSum[NVOXELS];
    for (int zz = 0; zz < NVOXELS; zz++) {
        temp[zz] = 0.f;
        if (no_norm == 0u)
            wSum[zz] = 0.f;
    }
    //const float nZ = 1.f / convert_float(d_nProjections);
    //float pz = nZ / 2.f;
    float3 dV = convert_float3(i) * d_d + d_d / 2.f + b;
    const float2 koko = { convert_float(d_size_x) * d_dPitch.x, convert_float(d_sizey) * d_dPitch.y };
    const float2 indeksi = { convert_float(d_size_x) / 2.f, convert_float(d_sizey) / 2.f };
    //float apuZ = d_d.z / 2.f + b.z;
//#pragma unroll NPROJECTIONS
    for (int kk = 0; kk < d_nProjections; kk++) {

        float3 d1, d2, d3;
        float3 s;
        //const float kulmaXY = d_angles[kk * NA];
        s = (float3)(d_xyz[kk * 6], d_xyz[kk * 6 + 1], d_xyz[kk * 6 + 2]);
        d1 = (float3)(d_xyz[kk * 6 + 3], d_xyz[kk * 6 + 4], d_xyz[kk * 6 + 5]);
        //s = vload3(kk * 2, d_xyz);
        //d1 = vload3(kk * 2, &d_xyz[3]);
#if defined(PITCH)
        const float3 apuX = (float3)(d_uv[kk * NA], d_uv[kk * NA + 1], d_uv[kk * NA + 2]) * indeksi.x;
        const float3 apuY = (float3)(d_uv[kk * NA + 3], d_uv[kk * NA + 4], d_uv[kk * NA + 5]) * indeksi.y;
        //d2.x = d1.x + indeksi.x * d_uv[kk * NA] + indeksi.y * d_uv[kk * NA + 1];
        //d3.x = d1.x - indeksi.x * d_uv[kk * NA] - indeksi.y * d_uv[kk * NA + 1];
        //d2.y = d1.y + indeksi.x * d_uv[kk * NA + 2] + indeksi.y * d_uv[kk * NA + 3];
        //d3.y = d1.y - indeksi.x * d_uv[kk * NA + 2] - indeksi.y * d_uv[kk * NA + 3];
        //d2.z = d1.z + indeksi.x * d_uv[kk * NA + 4] + indeksi.y * d_uv[kk * NA + 5];
        //d3.z = d1.z - indeksi.x * d_uv[kk * NA + 4] - indeksi.y * d_uv[kk * NA + 5];
#else
        const float3 apuX = (float3)(indeksi.x * d_uv[kk * NA], indeksi.x * d_uv[kk * NA + 1], 0.f);
        const float3 apuY = (float3)(0.f, 0.f, indeksi.y * d_dPitch.y);
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
        float lowerPart = dot(-v, cP);
        //const float3 dMin = d2 - d3;
        const float vApu = v.x * v.x + v.y * v.y;
        const float dApu = d_d.z * cP.z;
        const float pz = (convert_float(kk) + 0.5f) / convert_float(get_image_depth(d_forw));
#pragma unroll NVOXELS
        for (int zz = 0; zz < NVOXELS; zz++) {
            //d.z = dV.z;
            const uint ind = i.z + zz;
            if (ind >= d_N.z)
                break;
            const float t = native_divide(upperPart, lowerPart);
            //float3 p = s + v * t;
            float3 p = mad(v, t, s);
            const float L = distance(p, s);
            //const float l1 = v.x * v.x + v.y * v.y + v.z * v.z;
            const float l1 = vApu + v.z * v.z;
            p -= d3;
            const float weight = (L * L * L) / (l1)*kerroin;
            //const float yVar = read_imagef(d_forw, samplerIm, (float4)(dot(p, normX) / koko.x, dot(p, normY) / koko.y, pz, 0.f)).x;
            const float yVar = read_imagef(d_forw, samplerIm, (float4)(dot(p, normX), dot(p, normY), pz, 0.f)).x;
            //const float yVar = read_imagef(d_forw, sampler, (float4)(p.x, p.z, pz, 0.f)).x;
            //const float yVar = p.x;
            if (yVar != 0.f) {
                temp[zz] += yVar * weight;
                //temp += yVar;
                if (no_norm == 0u)
                    wSum[zz] += weight;
            }
            v.z += d_d.z;
            //lowerPart -= d_d.z * cP.z;
            lowerPart -= dApu;
        }
        //pz += nZ;
        //break;
//#endif
    }
#ifdef BP
    //if (idx < 428589460)
    for (int zz = 0; zz < NVOXELS; zz++) {
        const uint ind = i.z + zz;
        if (ind >= d_N.z)
            break;
        d_OSEM[idx] = temp[zz];
        if (no_norm == 0u)
            d_Summ[idx] = wSum[zz];
        idx += d_N.y * d_N.x;
    }
#endif
}
#endif