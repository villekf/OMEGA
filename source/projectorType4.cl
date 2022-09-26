
//__constant sampler_t sampler_MASK = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_NEAREST | CLK_ADDRESS_NONE;
#ifdef SPECTMASK
__constant sampler_t sampler_SPECT = CLK_NORMALIZED_COORDS_TRUE | CLK_FILTER_NEAREST | CLK_ADDRESS_CLAMP;
#endif

//#define NVOXELS 8
//#define PITCH

#if (defined(FP) && !defined(SPECTMASK))
#define NSTEPS 10000

// Forward projection
__kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, 1)))
void projectorType4Forward(const uint3 d_N, const float3 b, const uint d_size_x, const uint d_sizey, const float2 d_dPitch, const float3 bmax, 
    const float dL, const float3 d_scale,
#ifdef MASKFP
    __read_only image2d_t maskFP,
#endif
    __read_only image3d_t d_OSEM, __global float* d_forw, __constant float* d_xyz, __constant float* d_uv, const long d_nProjections)
{
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
#if !defined(CT) // PET and SPECT ONLY
#ifdef TOF // TOF ONLY
        float local_sino = 0.f;
#ifndef LISTMODE2
#pragma unroll NBINS
    for (long to = 0L; to < NBINS; to++)
        local_sino += d_Sino[idx + m_size * to];
#endif
#else
#ifdef LISTMODE2
        const float local_sino = 0.f;
#else
        const float local_sino = (d_Sino[idx]);
#endif
#endif
#ifndef MBSREM
        if (no_norm == 1u && local_sino == 0.f)
            return;
#else
        const uchar no_norm = 0u;
#endif
#endif

    float3 s, d;
#ifdef LISTMODE
    getDetectorCoordinatesListmode(d_xyz, &s, &d, idx);
#else
#if defined(CT) //|| (defined(SPECT) && !defined(SPECTMASK))
    getDetectorCoordinatesCT(d_xyz, d_uv, &s, &d, i, d_size_x, d_sizey, d_dPitch);
#elif defined(RAW) // raw data
    getDetectorCoordinatesRaw(d_xyz, d_zdet, d_L, d_det_per_ring, idx, &s, &d); // Sinogram data
#elif defined(FIND_LORS) || !defined(SUBSETS) // Precomputation phase
    getDetectorCoordinatesFullSinogram(d_size_x, i, &s, &d, d_xy, d_z);
#else // Not the precomputation phase
    getDetectorCoordinates(d_xyindex, d_zindex, d_size_x, idx, &s, &d, d_xy, d_zdet);
#endif
#endif
#if !defined(CT)
    bool RHS = false, SUMMA = false;
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
#ifdef MBSREM
#ifdef TOF
    float axACOSEM[NBINS];
    float ax[NBINS];
#pragma unroll NBINS
    for (uint to = 0; to < NBINS; to++) {
        axACOSEM[to] = 0.f;
        ax[to] = 0.f;
    }
#else
    float axACOSEM = 0.f;
    float axCOSEM = 0.f;
#endif
#ifdef TOF
    float minimi[NBINS];
#pragma unroll NBINS
    for (uint to = 0; to < NBINS; to++)
        minimi[to] = 1e8f;
#else
    float minimi = 1e8f;
#endif
#else
    float ax[NROLLS];
#pragma unroll
    for (uint kk = 0; kk < NROLLS; kk++)
        ax[kk] = 0.f;
#endif
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
    uint local_ind = 0u;
    float local_ele = 0.f;
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
    if (tStart >= tEnd)
        return;
#ifndef CT
    const float L = length(v);
    float jelppi = 0.f;
    const float tStep = dL / L;
    float dP = tStep;
#else
    //const float tStep = dL / length(v);
    const float tStep = native_divide(dL, length(v));
#endif

    //if (i.x == 1 && i.y == 300 && i.z == 59) {
    //    printf("s.x = %f\n", s.x);
    //    printf("s.y = %f\n", s.y);
    //    printf("s.z = %f\n", s.z);
    //    printf("d.x = %f\n", d.x);
    //    printf("d.y = %f\n", d.y);
    //    printf("d.z = %f\n", d.z);
    //    printf("bmin.x = %f\n", bmin.x);
    //    printf("bmin.y = %f\n", bmin.y);
    //    printf("bmin.z = %f\n", bmin.z);
    //    printf("bmax.x = %f\n", bmax.x);
    //    printf("bmax.y = %f\n", bmax.y);
    //    printf("bmax.z = %f\n", bmax.z);
    //    printf("NA = %d\n", NA);
    //    //printf("dPZ = %f\n", dPZ);
    //    //printf("dV.x = %f\n", dV.x);
    //    //printf("dV.y = %f\n", dV.y);
    //    //printf("dV.z = %f\n", dV.z);
    //}
    s = (s - bmin) * d_scale;
    v *= d_scale;
    //float3 p = (s + tStart * v - bmin) * d_scale;
    //float3 p = s + tStart * v;
#ifdef FP
    float temp = 0.f;

    float t = tStart;
    for (uint ii = 0; ii < NSTEPS; ii++) {
        const float4 p = (float4)(mad(t, v, s), 0.f);
#if !defined(CT) && !defined(AF)
#ifdef ATN
        jelppi -= read_imagef(d_atten, samplerIm, p).x;
#endif
#if defined(TOF) && (defined(DEC) || defined(FP))
        const float TOFSum = TOFLoop(DD, dL, store_elements, TOFCenter, sigma_x, &D, ii * NBINS, d_epps);
#endif
#ifdef FP
#if !defined(CT)
        if (local_sino != 0.f) {
#endif
#ifdef AF // Implementation 2

#ifdef MBSREM
#ifdef TOF
            if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0u)
                denominatorTOF(ax, dL, d_OSEM, p, TOFSum, store_elements, DD, TOFCenter, sigma_x, &D, ii * NBINS, d_epps, d_N);
#else
            if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && d_alku == 0u)
                denominator_multi(local_ele, &axCOSEM, read_imagef(d_OSEM, samplerIm, p).x);
#endif
#else
#ifdef TOF
            denominatorTOF(ax, dL, d_OSEM, p, TOFSum, store_elements, DD, TOFCenter, sigma_x, &D, ii * NBINS, d_epps, d_N);
#else
            denominator(dL, ax, p, d_N, d_OSEM);
#endif

#endif

#else // Implementation 3
#ifdef TOF
            denominatorTOF(ax, dL, d_OSEM, p, TOFSum, store_elements, DD, TOFCenter, sigma_x, &D, ii * NBINS, d_epps, d_N);
#else
            denominator_multi(dL, &axOSEM, read_imagef(d_OSEM, samplerIm, p).x);
#endif

#endif
#if !defined(CT)
        }
#endif
#endif
#else
        temp += read_imagef(d_OSEM, samplerIm, p).x;
#endif
        //if (i.x == 20 && i.y == 20 && (i.z) == 10 && ii == 0) {
        //    printf("p.y = %f\n", p.y);
        //    printf("p.x = %f\n", p.x);
        //    printf("p.z = %f\n", p.z);
        //}
        t += tStep;
        if (t > tEnd)
            break;
        //p = (s + t * v - bmin) * d_scale;
        //p = s + t * v;
        //p = mad(t, v, s);
    }
    //if (i.x == 20 && i.y == 20 && (i.z) == 10) {
    //    //printf("indeksi.x = %f\n", indeksi.x);
    //    //printf("indeksi.y = %f\n", indeksi.y);
    //    printf("i.x = %d\n", i.x);
    //    printf("i.y = %d\n", i.y);
    //    printf("i.z = %d\n", i.z);
    //    //printf("apuX.x = %f\n", d_uv[kk * NA]);
    //    //printf("apuX.y = %f\n", d_uv[kk * NA + 1]);
    //    //printf("apuX.z = %f\n", d_uv[kk * NA + 2]);
    //    //printf("apuY.x = %f\n", d_uv[kk * NA + 3]);
    //    //printf("apuY.y = %f\n", d_uv[kk * NA + 4]);
    //    //printf("apuY.z = %f\n", d_uv[kk * NA + 5]);
    //    //printf("bmin.x = %f\n", bmin.x);
    //    //printf("bmin.y = %f\n", bmin.y);
    //    //printf("bmin.z = %f\n", bmin.z);
    //    //printf("bmax.x = %f\n", bmax.x);
    //    //printf("bmax.y = %f\n", bmax.y);
    //    //printf("bmax.z = %f\n", bmax.z);
    //    printf("NA = %d\n", NA);
    //    printf("temp = %f\n", temp);
    //    printf("tStart = %f\n", tStart);
    //    printf("tEnd = %f\n", tEnd);
    //    //printf("zz = %d\n", zz);
    //    //printf("pz = %f\n", pz);
    //    //printf("d.z = %f\n", d.z);
    //    printf("s.x = %f\n", s.x);
    //    printf("s.y = %f\n", s.y);
    //    printf("s.z = %f\n", s.z);
    //    printf("d.x = %f\n", d.x);
    //    printf("d.y = %f\n", d.y);
    //    printf("d.z = %f\n", d.z);
    //}

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
#ifdef MBSREM
    if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM > 0) && local_sino != 0.f && d_alku == 0u) {
#ifdef TOF
#pragma unroll NBINS
        for (int to = 0; to < NBINS; to++) {
            ax[to] *= temp;
            if (ax[to] < d_epps)
                ax[to] = d_epps;
#ifdef RANDOMS
            ax[to] += d_sc_ra[idx];
#endif
            ax[to] = d_Sino[idx + to * m_size] / ax[to];
        }
#else
#ifndef CT
        axCOSEM *= temp;
        if (axCOSEM < d_epps)
            axCOSEM = d_epps;
#endif
#ifdef RANDOMS
        axCOSEM += d_sc_ra[idx];
#endif
        axCOSEM = local_sino / axCOSEM;
#endif
    }
    RHS = true;
#else
#ifndef AF
#ifdef FP 
    if (fp == 1) {
#ifdef TOF
        nominatorTOF(MethodList, ax, d_Sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx, m_size, local_sino);
#pragma unroll NBINS
        for (int to = 0; to < NBINS; to++)
            d_forw[idx + to * m_size] = ax[to];
#else
        nominator_multi(&axOSEM, local_sino, d_epps, temp, d_sc_ra, idx);
        d_forw[idx] = axOSEM;
#endif
        return;
    }
#endif
#endif

#ifdef FP
    if (local_sino != 0.f) {
#ifdef TOF
        nominatorTOF(MethodList, ax, d_Sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx, m_size, local_sino);
#else
#ifdef AF
        nominator(MethodList, ax, local_sino, d_epsilon_mramla, d_epps, temp, d_sc_ra, idx);
#else
        nominator_multi(&axOSEM, local_sino, d_epps, temp, d_sc_ra, idx);
#endif
#endif
    }
#endif
    RHS = true;
    //else
    //	SUMMA = true;
#endif
    t = tStart;
    const float4 Nf = (float4)(convert_float3(d_N - 1), 0.f);
    for (uint ii = 0; ii < NSTEPS; ii++) {
        const float4 p = (float4)(mad(t, v, s), 0.f);
        const int4 ind = convert_int4_rtz((p * Nf));
        const uint local_ind = ind.x + ind.y * d_N.x + ind.z * d_N.x * d_N.y;
#ifdef TOF
#ifndef DEC
        const float TOFSum = TOFLoop(DD, dL / temp, store_elements, TOFCenter, sigma_x, &D, ii * NBINS, d_epps);
#endif
        backprojectTOF(local_ind, p, dL, ii * NBINS, store_elements, ax, d_Summ, local_sino,
#ifndef DEC
            temp, sigma_x, &D, DD, TOFCenter, d_epps, TOFSum,
#endif
#ifdef MBSREM
            MethodListOpenCL, d_alku, MBSREM_prepass, minimi, axACOSEM, d_OSEM, d_E, d_co, d_aco, idx, m_size);
#else
            d_forw, no_norm, d_N);
#endif
#else
#ifdef MBSREM
        if (d_alku == 0u) {
            if (MBSREM_prepass == 1)
#ifdef ATOMIC
                atom_add(&d_Summ[local_ind], convert_long(dL * TH));
#elif defined(ATOMIC32)
                atomic_add(&d_Summ[local_ind], convert_int(dL * TH));
#else
                atomicAdd_g_f(&d_Summ[local_ind], dL);
#endif
            if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1) {
                if (local_ele < minimi && local_ele > 0.f)
                    minimi = local_ele;
                d_E[idx] += local_ele;
            }
            if ((MethodListOpenCL.COSEM == 1 || MethodListOpenCL.ECOSEM == 1 || MethodListOpenCL.OSLCOSEM == 2) && local_sino != 0.f)
#ifdef ATOMIC
                atom_add(&d_co[local_ind], convert_long(axCOSEM * dL * TH));
#elif defined(ATOMIC32)
                atomic_add(&d_co[local_ind], convert_int(axCOSEM * dL * TH));
#else
                atomicAdd_g_f(&d_co[local_ind], axCOSEM * dL);
#endif
            if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && local_sino != 0.f)
#ifdef ATOMIC
                atom_add(&d_aco[local_ind], convert_long(axCOSEM * dL * TH));
#elif defined(ATOMIC32)
                atomic_add(&d_aco[local_ind], convert_int(axCOSEM * dL * TH));
#else
                atomicAdd_g_f(&d_aco[local_ind], axCOSEM * dL);
#endif
        }
        if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u)
            axACOSEM += (dL * read_imagef(d_OSEM, samplerIm, p).x);
        //axACOSEM += (local_ele * d_OSEM[local_ind]);
#else
        if (no_norm == 0u)
#ifdef ATOMIC
            atom_add(&d_Summ[local_ind], convert_long(dL * TH));
#elif defined(ATOMIC32)
            atomic_add(&d_Summ[local_ind], convert_int(dL * TH));
#else
            atomicAdd_g_f(&d_Summ[local_ind], dL);
#endif

        if (local_sino != 0.f) {
#ifdef AF
            rhs(MethodList, dL, ax, local_ind, d_N, d_forw);
#else

#ifdef ATOMIC
            atom_add(&d_forw[local_ind], convert_long(dL * axOSEM * TH));
#elif defined(ATOMIC32)
            atomic_add(&d_forw[local_ind], convert_int(dL * axOSEM * TH));
#else
            atomicAdd_g_f(&d_forw[local_ind], (dL * axOSEM));
#endif
#endif
        }
#endif
#endif
    }
#ifdef MBSREM
#ifdef TOF
#pragma unroll NBINS
    for (long to = 0L; to < NBINS; to++) {
        if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1 && d_alku == 0u)
            d_Amin[idx + to * m_size] = minimi[to];
        if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u) {
#ifdef RANDOMS
            axACOSEM[to] += d_sc_ra[idx];
#endif
            d_ACOSEM_lhs[idx + to * m_size] = axACOSEM[to];
        }
    }
#else
    if ((MethodListOpenCL.MRAMLA_ == 1 || MethodListOpenCL.MBSREM_ == 1) && MBSREM_prepass == 1 && d_alku == 0u)
        d_Amin[idx] = minimi;
    if ((MethodListOpenCL.ACOSEM == 1 || MethodListOpenCL.OSLCOSEM == 1) && d_alku > 0u) {
#ifdef RANDOMS
        axACOSEM += d_sc_ra[idx];
#endif
        d_ACOSEM_lhs[idx] = axACOSEM;
    }
#endif
#endif
#else
    d_forw[idx] = temp * dL;
#endif
#endif

}

#endif

#if defined(BP) || (defined(SPECTMASK) && defined(FP))
//#define MAX(a,b) (a>b?a:b)
//#define MIN(a,b) (a<b?a:b)
//#undef PITCH
//#undef NA
//#define NA 6
__kernel __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, 1)))
void projectorType4Backward(const uint3 d_N, const float3 b, const uint d_size_x, const uint d_sizey, const float2 d_dPitch, const float3 d_d, const float kerroin,
#ifdef SPECTMASK
    const float cThickness, __read_only image2d_t d_mask1, __read_only image2d_t d_mask2,
#endif
#ifdef MASKBP
    __read_only image2d_t maskBP,
#endif
    __read_only image3d_t d_forw, __global float* d_OSEM, __global float* d_Summ, __constant float* d_xyz, __constant float* d_uv,
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
    float3 dV = convert_float(i) * d_d + d_d / 2.f + b;
    const float2 koko = { convert_float(d_size_x) * d_dPitch.x, convert_float(d_sizey) * d_dPitch.y };
    //float apuZ = d_d.z / 2.f + b.z;
//#pragma unroll NPROJECTIONS
    for (int kk = 0; kk < d_nProjections; kk++) {

        float3 d1, d2, d3;
        float3 s;
        //const float kulmaXY = d_angles[kk * NA];
#ifdef SPECTMASK
            s = d;
            dd.x = d_x[kk];
            dd.y = d_y[kk];
            dd.z = d_zdet[kk];
            float3 d1, d2, d3, apu, v, p, dD;
            float t;
            //int ii = 0;
            //int jj = 0;
            for (int ii = 0; ii < size_x; ii++) {
                dD.x = dd.x - d_dPitch * convert_float(ii) * native_cos(kulmaXY);
                dD.y = dd.y - d_dPitch * convert_float(ii) * native_sin(kulmaXY);
                for (int jj = 0; jj < size_z; jj++) {
                    dD.z = dd.z + d_dPitch * convert_float(jj);
                    v = dD - s;
#if defined(MASK2)
                    d1 = (float3)(dd.x - cThickness * native_sin(kulmaXY), dd.y - cThickness * native_cos(kulmaXY), dd.z);
                    d2 = (float3)(d1.x, d1.y, d1.z + d_dPitch * 10.f);
                    d3 = (float3)(d1.x - d_dPitch * native_cos(kulmaXY) * 10.f, d1.y - d_dPitch * native_sin(kulmaXY) * 10.f, d1.z);
                    apu = cross(d2 - d1, d3 - d1);
                    t = dot(apu, s - d1) / dot(-v, apu);
                    p = s + v * t;
                    p.x = ((d1.x * native_cos(M_PI_F * 2.f - kulmaXY) - d1.y * native_sin(M_PI_F * 2.f - kulmaXY)) - (p.x * native_cos(M_PI_F * 2.f - kulmaXY) - p.y * native_sin(M_PI_F * 2.f - kulmaXY))) / (d_dSize.y);
                    p.z -= (d1.z - d_dPitch / 2.f);
                    p.z /= d_dSize.z;
                    const float maskVal1 = read_imagef(d_mask1, sampler_SPECT, (float2)(p.x, p.z)).x;
                    //temp = maskVal1;
                    if (maskVal1 == 0.f)
                        continue;
#endif
                    //d1 = (float3)(dd.x - d_dPitch * native_cos(kulmaXY) * 10.f, dd.y - d_dPitch * native_sin(kulmaXY) * 10.f, dd.z );
                    //d2 = (float3)(dd.x, dd.y, dd.z + d_dPitch * 10.f );
                    //apu = cross(d1 - dd, d2 - dd);
                    //t = dot(apu, s - dd) / dot(-v, apu);
                    //p = s + v * t;
                    p = dD;
                    const float L = distance(p, s);
                    p.x = ((dd.x * native_cos(M_PI_F * 2.f - kulmaXY) - dd.y * native_sin(M_PI_F * 2.f - kulmaXY)) - (p.x * native_cos(M_PI_F * 2.f - kulmaXY) - p.y * native_sin(M_PI_F * 2.f - kulmaXY))) / (d_dSize.y);
                    //p.x = (p.x * native_cos(M_PI_F * 2.f - kulmaXY) - p.y * native_sin(M_PI_F * 2.f - kulmaXY)) / (d_dSize.y);
                    p.z -= (dd.z - d_dPitch / 2.f);
                    //p.z += d_dPitch / 2.f;
                    p.z /= d_dSize.z;
                    //const float maskVal2 = read_imagef(d_mask2, sampler_SPECT, (float2)(p.x, p.z)).x;
                    //if (maskVal2 == 0.f)
                    //    continue;
                    //const float l = v.x * v.x + v.y * v.y + v.z * v.z;
                    //const float weight = (L * L * L) / (l)*kerroin;
                    const float weight = 1.f;
#ifdef BP
                    const float yVar = read_imagef(d_forw, sampler_SPECT, (float4)(p.x, p.z, pz + nZ * kk, 0.f)).x;
#else
                    const float yVar = read_imagef(d_forw, sampler_SPECT, (float4)((d.x - b.x) / (convert_float(d_Nx) * d_dx), (d.y - b.y) / (convert_float(d_Ny) * d_dy), (d.z - b.z) / (convert_float(d_Nz) * d_dz), 0.f)).x;
                    //const float yVar = maskVal1;
#endif
                //temp = p.x;
                    if (yVar > 0.f) {
#ifdef BP
                        temp += yVar * (weight / L);
                        if (no_norm == 0u)
                            wSum += (weight / L);
#else
                        atomicAdd_g_f(&d_OSEM[kk * size_z * size_x + ii + jj * size_x], yVar * (weight / L));
#endif
                    }
                }
            }
#else
        s = (float3)(d_xyz[kk * 6], d_xyz[kk * 6 + 1], d_xyz[kk * 6 + 2]);
        d1 = (float3)(d_xyz[kk * 6 + 3], d_xyz[kk * 6 + 4], d_xyz[kk * 6 + 5]);
        //s = vload3(kk * 2, d_xyz);
        //d1 = vload3(kk * 2, &d_xyz[3]);
        const float2 indeksi = { convert_float(d_size_x) / 2.f, convert_float(d_sizey) / 2.f };
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
        const float3 normX = normalize(apuX);
        const float3 normY = normalize(apuY);
        const float3 cP = cross(d2, d3 - d1);
        const float upperPart = dot(cP, s - d1);
        float3 v = dV - s;
        float lowerPart = dot(-v, cP);
        //const float3 dMin = d2 - d3;
#endif
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
            const float l1 = v.x * v.x + v.y * v.y + v.z * v.z;
            p -= d3;
            const float weight = (L * L * L) / (l1)*kerroin;
            const float yVar = read_imagef(d_forw, samplerIm, (float4)(dot(p, normX) / koko.x, dot(p, normY) / koko.y, pz, 0.f)).x;
            //const float yVar = read_imagef(d_forw, sampler, (float4)(p.x, p.z, pz, 0.f)).x;
            //if (i.x == 20 && i.y == 20 && (ind) == 10 && kk == 20) {
            //    printf("indeksi.x = %f\n", indeksi.x);
            //    printf("indeksi.y = %f\n", indeksi.y);
            //    printf("i.x = %d\n", i.x);
            //    printf("i.y = %d\n", i.y);
            //    printf("i.z = %d\n", i.z);
            //    //printf("p.y = %f\n", p.y);
            //    //printf("p.x = %f\n", p.x);
            //    //printf("p.z = %f\n", p.z);
            //    printf("apuX.x = %f\n", d_uv[kk * NA]);
            //    printf("apuX.y = %f\n", d_uv[kk * NA + 1]);
            //    printf("apuX.z = %f\n", d_uv[kk * NA + 2]);
            //    printf("apuY.x = %f\n", d_uv[kk * NA + 3]);
            //    printf("apuY.y = %f\n", d_uv[kk * NA + 4]);
            //    printf("apuY.z = %f\n", d_uv[kk * NA + 5]);
            //    //printf("v.x = %f\n", v.x);
            //    //printf("v.y = %f\n", v.y);
            //    //printf("v.z = %f\n", v.z);
            //    //printf("NA = %d\n", NA);
            //    printf("yVar = %f\n", yVar);
            //    printf("zz = %d\n", zz);
            //    //printf("pz = %f\n", pz);
            //    //printf("d.z = %f\n", d.z);
            //    printf("s.x = %f\n", s.x);
            //    printf("s.y = %f\n", s.y);
            //    printf("s.z = %f\n", s.z);
            //    printf("d.x = %f\n", d1.x);
            //    printf("d.y = %f\n", d1.y);
            //    printf("d.z = %f\n", d1.z);
            //    //printf("dPZ = %f\n", dPZ);
            //    printf("dot(p, normX) / koko.x = %f\n", dot(p, normX) / koko.x);
            //    printf("dot(p, normX) / koko.x = %f\n", dot(p, normY) / koko.y);
            //    //printf("dV.y = %f\n", dV.y);
            //    //printf("dV.z = %f\n", dV.z);
            //}
            //const float yVar = p.x;
            if (yVar != 0.f) {
                temp[zz] += yVar * weight;
                //temp += yVar;
                if (no_norm == 0u)
                    wSum[zz] += weight;
            }
            v.z += d_d.z;
            lowerPart -= d_d.z * cP.z;
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


//#ifndef AF
//__kernel void summa(const __global CAST* d_Summ_device, __global CAST* d_Summ_local, const __global CAST* d_rhs_device, __global CAST* d_rhs_local,
//    const uint im_dim, const uchar no_norm) {
//
//    uint gid = get_global_id(0);
//
//    for (uint i = gid; i < im_dim; i += get_global_size(0)) {
//        if (no_norm == 0u)
//            d_Summ_local[i] += d_Summ_device[i];
//        d_rhs_local[i] += d_rhs_device[i];
//    }
//}
//#endif