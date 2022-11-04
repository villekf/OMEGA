
__constant sampler_t sampler2 = CLK_NORMALIZED_COORDS_TRUE | CLK_FILTER_LINEAR | CLK_ADDRESS_CLAMP_TO_EDGE;

//#if (defined(FP) && !defined(SPECTMASK))
//#define NSTEPS 10000

// Forward projection for the branchless distance-driven ray tracer
#ifdef FP
//#undef PITCH
//#undef NA
//#define NA 6
__kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, 1)))
void projectorType5Forward(const uint3 d_N, const float3 b, const uint d_size_x, const uint d_sizey, const float2 d_dPitch, const float3 d_d,
    const float3 d_scale, const float3 d_Size,
#ifdef MASKFP
    __read_only image2d_t maskFP,
#endif
    __constant float* d_xyz, __constant float* d_uv, __read_only image3d_t d_IImageY, __read_only image3d_t d_IImageX, __global float* d_forw,
#ifdef MEANDISTANCEFP
    __constant float* d_meanV,
#endif
    const long d_nProjections) {

    const int3 i = { get_global_id(0), get_global_id(1), get_global_id(2) };
    size_t idx = get_global_id(0) + get_global_id(1) * d_size_x + get_global_id(2) * d_sizey * d_size_x;

    if (i.x >= d_size_x || i.y >= d_sizey || i.z >= d_nProjections)
        return;

#ifdef MASKFP
    const int maskVal = read_imageui(maskFP, sampler_MASK, (int2)(i.x, i.y)).x;
    if (maskVal == 0)
        return;
#endif
    float3 s, d, dL, dR, dU, dD;
    //float3 s, d;
    float temp = 0.f;
#ifdef CT
    // Get source and detector coordinates
    getDetectorCoordinatesCT(d_xyz, d_uv, &s, &d, i, d_size_x, d_sizey, d_dPitch, &dL, &dR, &dU, &dD);
    // Weight values
    const float3 uVector = fabs(normalize(d - s));
    float kerroin = d_d.x * d_d.y * d_d.z;
    const float apuZ = (b.z - s.z - d_d.z / 2.f);
    //const float apuZ = (b.z - s.z - d_d.z);
    if (fabs(s.x) <= fabs(s.y)) {
        const float2 X = { (dL.x - s.x) / (dL.y - s.y), (dR.x - s.x) / (dR.y - s.y) };
        const float2 Z = { (dU.z - s.z) / (dU.y - s.y), (dD.z - s.z) / (dD.y - s.y) };
        kerroin /= uVector.y;
        //const float apuX = (b.x - s.x - d_d.x);
        const float apuX = (b.x - s.x - d_d.x / 2.f);
        const float apuV = b.y - s.y;
        for (uint jj = 0; jj < d_N.y; jj++) {
            float dy = convert_float(jj) * d_d.y + d_d.y / 2.f;
            const float v = dy + apuV;
            float2 xLR = v * X;
            float2 zUD = v * Z;
            const float area = fabs((xLR.x - xLR.y) * (zUD.x - zUD.y));
            xLR -= apuX;
            zUD -= apuZ;
            xLR *= d_scale.x;
            dy *= d_Size.y;
            zUD *= d_scale.z;
//#ifdef MEANDISTANCEFP
//            if (any(xLR > 1.f) || any(zUD > 1.f) || any(xLR < 0.f) || any(zUD < 0.f))
//                continue;
//#endif
            if (xLR.x > xLR.y) {
                const float apu = xLR.x;
                xLR.x = xLR.y;
                xLR.y = apu;
            }
            if (zUD.x > zUD.y) {
                const float apu = zUD.x;
                zUD.x = zUD.y;
                zUD.y = apu;
            }
            const float A = read_imagef(d_IImageY, sampler2, (float4)(xLR.y, zUD.x, dy, 0.f)).x;
            const float B = read_imagef(d_IImageY, sampler2, (float4)(xLR.x, zUD.y, dy, 0.f)).x;
            const float C = read_imagef(d_IImageY, sampler2, (float4)(xLR.y, zUD.y, dy, 0.f)).x;
            const float D = read_imagef(d_IImageY, sampler2, (float4)(xLR.x, zUD.x, dy, 0.f)).x;
            float apu = C + D - A - B;
            //if (i.x == 300 && i.y == 300 && i.z == 0 && jj < 10) {
            //    printf("apu1 = %f\n", apu);
            //}
            //if (i.x == 300 && i.y == 300 && i.z == 0 && jj == 10) {
            //    printf("dL.x = %f\n", dL.x);
            //    printf("dL.y = %f\n", dL.y);
            //    printf("dR.x = %f\n", dR.x);
            //    printf("dR.y = %f\n", dR.y);
            //    printf("dU.y = %f\n", dU.y);
            //    printf("dU.z = %f\n", dU.z);
            //    printf("dD.y = %f\n", dD.y);
            //    printf("dD.z = %f\n", dD.z);
            //    printf("apu1 = %f\n", apu);
            //    printf("C = %f\n", C);
            //    printf("D = %f\n", D);
            //    printf("A = %f\n", A);
            //    printf("B = %f\n", B);
            //}
#ifdef MEANDISTANCEFP
            const float area2 = fabs((xLR.x - xLR.y) * (zUD.x - zUD.y)) * convert_float(get_image_width(d_IImageY) * get_image_height(d_IImageY));
            apu += d_meanV[jj + d_N.x] * area2;
#endif
            temp += apu / area;
            //if (i.x == 300 && i.y == 1 && i.z == 49) {
            //    printf("apu = %f\n", apu);
            //    printf("d_meanV[jj + d_N.x] = %f\n", d_meanV[jj + d_N.x]);
            //    printf("area2 = %f\n", area2);
            //}
        }
    }
    else {
        const float2 Y = { (dL.y - s.y) / (dL.x - s.x), (dR.y - s.y) / (dR.x - s.x) };
        const float2 Z = { (dU.z - s.z) / (dU.x - s.x), (dD.z - s.z) / (dD.x - s.x) };
        kerroin /= uVector.x;
        //const float apuY = (b.y - s.y - d_d.y);
        const float apuY = (b.y - s.y - d_d.y / 2.f);
        const float apuV = b.x - s.x;
        for (uint ii = 0; ii < d_N.x; ii++) {
            float dx = convert_float(ii) * d_d.x + d_d.x / 2.f;
            const float v = dx + apuV;
            float2 yLR = v * Y;
            float2 zUD = v * Z;
            const float area = fabs((yLR.x - yLR.y) * (zUD.x - zUD.y));
            yLR -= apuY;
            zUD -= apuZ;
            yLR *= d_scale.y;
            dx *= d_Size.x;
            zUD *= d_scale.z;
//#ifdef MEANDISTANCEFP
//            if (any(yLR > 1.f) || any(zUD > 1.f) || any(yLR < 0.f) || any(zUD < 0.f))
//                continue;
//#endif
            if (yLR.x > yLR.y) {
                const float apu = yLR.x;
                yLR.x = yLR.y;
                yLR.y = apu;
            }
            if (zUD.x > zUD.y) {
                const float apu = zUD.x;
                zUD.x = zUD.y;
                zUD.y = apu;
            }
            const float A = read_imagef(d_IImageX, sampler2, (float4)(yLR.y, zUD.x, dx, 0.f)).x;
            const float B = read_imagef(d_IImageX, sampler2, (float4)(yLR.x, zUD.y, dx, 0.f)).x;
            const float C = read_imagef(d_IImageX, sampler2, (float4)(yLR.y, zUD.y, dx, 0.f)).x;
            const float D = read_imagef(d_IImageX, sampler2, (float4)(yLR.x, zUD.x, dx, 0.f)).x;
            float apu = C + D - A - B;
            //float apuT = apu;
            //if (i.x == 1 && i.y == 300 && i.z == 149) {
            //    printf("apu1 = %f\n", apu);
            //}
#ifdef MEANDISTANCEFP
            const float area2 = fabs((yLR.x - yLR.y) * (zUD.x - zUD.y)) * convert_float(get_image_width(d_IImageX) * get_image_height(d_IImageX));
            //const float testi = apu + d_meanV[ii] * area * convert_float(get_image_width(d_IImageX) * get_image_height(d_IImageX)) * d_scale.y * d_scale.z;
            //const float areaT = area * d_scale.y * d_scale.z;
            apu += (d_meanV[ii] * area2);
#endif
            temp += apu / area;
            //temp += (C + D - A - B) / area;
            //if (i.x == 1 && i.y == 300 && i.z == 149 && ii < 300 && ii > 250) {
            //    printf("temp = %f\n", temp);
            //    printf("apu = %f\n", apu);
            //    printf("apuT = %f\n", apuT);
            //    printf("testi = %f\n", testi);
            //    printf("d_meanV[ii] = %f\n", d_meanV[ii]);
            //    printf("area = %.8f\n", area);
            //    printf("area2 = %.8f\n", area2);
            //    printf("areaT = %.8f\n", areaT);
            //    printf("yLR.x = %f\n", yLR.x);
            //    printf("yLR.y = %f\n", yLR.y);
            //    printf("zUD.x = %f\n", zUD.x);
            //    printf("zUD.y = %f\n", zUD.y);
            //    printf("ii = %d\n", ii);
            //}
        }
    }
    //if (i.x == 1 && i.y == 300 && i.z == 149) {
    //if (i.x == 300 && i.y == 300 && i.z == 0) {
    //    printf("temp * kerroin = %f\n", temp * kerroin);
    //    printf("temp = %f\n", temp);
    //    printf("kerroin = %f\n", kerroin);
    //}
    d_forw[idx] = temp * kerroin;
#else

#endif
}
#endif

#ifdef BP

//#undef PITCH
//#undef NA
//#define NA 6
__kernel __attribute__((vec_type_hint(float))) __attribute__((reqd_work_group_size(LOCAL_SIZE, LOCAL_SIZE2, 1)))
void projectorType5Backward(const uint3 d_N, const float3 b, const uint d_size_x, const uint d_sizey, const float2 d_dPitch, const float3 d_d,
    const float3 d_scale, const float2 d_Size, 
#ifdef MASKBP
    __read_only image2d_t maskBP,
#endif
    __constant float* d_xyz, __constant float* d_uv, __read_only image3d_t d_IImage, __global float* d_forw, __global float* d_Summ,
#ifdef MEANDISTANCEBP
    __constant float* d_meanV,
#endif
    const uchar no_norm, const long d_nProjections) {

    const uint3 i = { get_global_id(0), get_global_id(1), get_global_id(2) };
    const size_t idx = get_global_id(0) + get_global_id(1) * d_N.x + get_global_id(2) * d_N.y * d_N.x;

    if (i.x >= d_N.x || i.y >= d_N.y || i.z >= d_N.z)
        return;

#ifdef MASKBP
    const int maskVal = read_imageui(maskBP, sampler_MASK, (int2)(i.x, i.y)).x;
    if (maskVal == 0)
        return;
#endif
    //float3 s, d;
    //float3 d2 = { 0.f, 0.f, 0.f };
    //float3 d3 = { 0.f, 0.f, 0.f };
    float temp = 0.f;
    float wSum = 0.f;
    //const float nZ = 1.f / convert_float(d_nProjections);
    //float pz = nZ / 2.f;
    const float3 dV = convert_float(i) * d_d + d_d / 2.f + b;
    const float2 koko = { convert_float(get_image_width(d_IImage)) * d_dPitch.x, convert_float(get_image_height(d_IImage)) * d_dPitch.y };
    const float2 indeksi = { convert_float(get_image_width(d_IImage)) / 2.f, convert_float(get_image_height(d_IImage)) / 2.f };
#ifdef CT
    for (int kk = 0; kk < d_nProjections; kk++) {

        float3 vLU, vLD, vRU, vRD;
        float3 d, d2, d3;
        //float3 s;
        int id = kk * 6;
        const float3 s = (float3)(d_xyz[id], d_xyz[id + 1], d_xyz[id + 2]);
        d = (float3)(d_xyz[id + 3], d_xyz[id + 4], d_xyz[id + 5]);
        //const float2 indeksi = { ceil(convert_float(get_image_width(d_IImage)) / 2.f + .4f), ceil(convert_float(get_image_height(d_IImage)) / 2.f + .4f) };
//#if defined(PITCH)
//        const float3 apuX = (float3)(d_uv[kk * NA], d_uv[kk * NA + 1], d_uv[kk * NA + 2]) * indeksi.x;
//        const float3 apuY = (float3)(d_uv[kk * NA + 3], d_uv[kk * NA + 4], d_uv[kk * NA + 5]) * indeksi.y;
//#else
        const float3 apuX = (float3)(indeksi.x * d_uv[kk * NA], indeksi.x * d_uv[kk * NA + 1], 0.f);
        const float3 apuY = (float3)(0.f, 0.f, indeksi.y * d_dPitch.y);
//#endif
        d2 = apuX - apuY;
        d3 = d - apuX - apuY;
        const float3 normX = normalize(apuX);
        const float3 normY = normalize(apuY);
        //const float2 v = dV.xy - s.xy;
        const float3 v = normalize(dV - s);
        //float3 v = (dV - s);
        const float3 crossP = cross(d2, d3 - d);
        const float upperPart = dot(crossP, s - d);
        if (fabs(s.x) <= fabs(s.y)) {
            if (s.y > 0.f) {
                vLU = (float3)(dV.x + d_d.x / 2.f, dV.y, dV.z + d_d.z / 2.f) - s;
                vLD = (float3)(dV.x + d_d.x / 2.f, dV.y, dV.z - d_d.z / 2.f) - s;
                vRU = (float3)(dV.x - d_d.x / 2.f, dV.y, dV.z + d_d.z / 2.f) - s;
                vRD = (float3)(dV.x - d_d.x / 2.f, dV.y, dV.z - d_d.z / 2.f) - s;
            }
            else {
                vLD = (float3)(dV.x + d_d.x / 2.f, dV.y, dV.z + d_d.z / 2.f) - s;
                vLU = (float3)(dV.x + d_d.x / 2.f, dV.y, dV.z - d_d.z / 2.f) - s;
                vRD = (float3)(dV.x - d_d.x / 2.f, dV.y, dV.z + d_d.z / 2.f) - s;
                vRU = (float3)(dV.x - d_d.x / 2.f, dV.y, dV.z - d_d.z / 2.f) - s;
            }
        }
        else {
            if (s.x > 0.f) {
                vLD = (float3)(dV.x, dV.y + d_d.y / 2.f, dV.z + d_d.z / 2.f) - s;
                vLU = (float3)(dV.x, dV.y + d_d.y / 2.f, dV.z - d_d.z / 2.f) - s;
                vRD = (float3)(dV.x, dV.y - d_d.y / 2.f, dV.z + d_d.z / 2.f) - s;
                vRU = (float3)(dV.x, dV.y - d_d.y / 2.f, dV.z - d_d.z / 2.f) - s;
            }
            else {
                vLU = (float3)(dV.x, dV.y + d_d.y / 2.f, dV.z + d_d.z / 2.f) - s;
                vLD = (float3)(dV.x, dV.y + d_d.y / 2.f, dV.z - d_d.z / 2.f) - s;
                vRU = (float3)(dV.x, dV.y - d_d.y / 2.f, dV.z + d_d.z / 2.f) - s;
                vRD = (float3)(dV.x, dV.y - d_d.y / 2.f, dV.z - d_d.z / 2.f) - s;
            }
            //continue;
        }
//#ifdef PITCH
//        const float tLU = upperPart / (dot(-vLU, crossP));
//        const float tRU = upperPart / (dot(-vRU, crossP));
//        const float tLD = upperPart / (dot(-vLD, crossP));
//        const float tRD = upperPart / (dot(-vRD, crossP));
//        float3 pLU = mad(vLU, tLU, s);
//        float3 pRU = mad(vRU, tRU, s);
//        float3 pLD = mad(vLD, tLD, s);
//        float3 pRD = mad(vRD, tRD, s);
//#else
        const float t = upperPart / (dot(-vLU, crossP));
        float3 pLU = mad(vLU, t, s);
        float3 pRU = mad(vRU, t, s);
        float3 pLD = mad(vLD, t, s);
        float3 pRD = mad(vRD, t, s);
//#endif
        pLU -= d3;
        pRU -= d3;
        pLD -= d3;
        pRD -= d3;
        const float dz = (convert_float(kk) + 0.5f) / convert_float(get_image_depth(d_IImage));
        const float4 coordA = { dot(pLU, normX) / koko.x, dot(pLU, normY) / koko.y, dz, 0.f };
        const float4 coordB = { dot(pRU, normX) / koko.x, dot(pRU, normY) / koko.y, dz, 0.f };
        const float4 coordC = { dot(pLD, normX) / koko.x, dot(pLD, normY) / koko.y, dz, 0.f };
        const float4 coordD = { dot(pRD, normX) / koko.x, dot(pRD, normY) / koko.y, dz, 0.f };
        //#ifdef MEANDISTANCEBP
        //if (any(coordA > 1.f) || any(coordB > 1.f) || any(coordA < 0.f) || any(coordB < 0.f) || any(coordC > 1.f) || any(coordD > 1.f) || any(coordC < 0.f) || any(coordD < 0.f))
        //    continue;
        //#endif
        float A = read_imagef(d_IImage, sampler2, coordA).x;
        float B = read_imagef(d_IImage, sampler2, coordB).x;
        float C = read_imagef(d_IImage, sampler2, coordC).x;
        float D = read_imagef(d_IImage, sampler2, coordD).x;
        //const float l = native_sqrt(v.x * v.x + v.y * v.y);
        //const float l = length(v);
        //v = normalize(v);
        float kerroin1;
        const float area = fabs(coordA.x - coordD.x) * fabs(coordA.y - coordD.y) * convert_float(get_image_width(d_IImage) * get_image_height(d_IImage));
        if (fabs(s.x) <= fabs(s.y))
            kerroin1 = fabs(v.y);
        else
            kerroin1 = fabs(v.x);
        const float kerroin = d_d.x / (kerroin1);
        //if (i.z == 100 && i.x == 200 && i.y == 10) {
        //    //printf("l = %f", l);
        //    printf("kerroin1 = % f", kerroin1);
        //    printf("  wSum = %f", wSum);
        //    printf("  kerroin = %f\n", kerroin);
        //    //printf("(d_d.x * d_d.y * d_d.z) = %f", (d_d.x * d_d.y * d_d.z));
        //    //printf("  v.y = %f", v.y);
        //    printf("area = %f\n", area);
        //    //printf("  v.z = %f\n", v.z);
        //}
        //float jelppi = A + D - C - B;
        //if (D < C) {
        //    const float ap = D;
        //    D = C;
        //    C = ap;
        //}
        //if (A > B) {
        //    const float ap = A;
        //    A = B;
        //    B = ap;
        //}
        //if (D > C && D < 0.f) {
        //    const float ap = D;
        //    D = C;
        //    C = ap;
        //}
        //if (A < B && A < 0.f) {
        //    const float ap = A;
        //    A = B;
        //    B = ap;
        //}
        float apu = A + D - C - B;
        //if (((B < A && B < D) || (C < A && C < D) || (C > A && C > D) || (B > A && B > D))) {
        ////if (((B < A && B < D) || (C < A && C < D) || (C > A && C > D) || (B > A && B > D)) && i.z == 8 && i.x == 165 && i.y == 15) {
        ////if (i.z == 79 && i.x == 127 && i.y == 35 && jelppi < 0.f) {
        ////if (i.z == 79 && i.x == 100 && i.y == 120) {
        //    	printf("A = %f", A);
        //    	printf("  B = %f", B);
        //    	printf("  C = %f", C);
        //    	printf("  D = %f\n", D);
        //    	//printf("  jelppi = %f", jelppi);
        //    	//printf("  apu = %f\n", apu);
        //    	//printf("coordA.x = %f", coordA.x * convert_float(get_image_width(d_IImage)));
        //    	//printf(" coordA.y = %f", coordA.y * convert_float(get_image_height(d_IImage)));
        //    	//printf(" coordA.z = %f\n", coordA.z);
        //    	//printf("coordB.x = %f", coordB.x* convert_float(get_image_width(d_IImage)));
        //    	//printf(" coordB.y = %f", coordB.y* convert_float(get_image_height(d_IImage)));
        //    	//printf(" coordB.z = %f\n", coordB.z);
        //    	//printf("coordC.x = %f", coordC.x);
        //    	//printf(" coordC.y = %f", coordC.y);
        //    	//printf(" coordC.z = %f\n", coordC.z);
        //    	//printf("coordD.x = %f", coordD.x);
        //    	//printf(" coordD.y = %f", coordD.y);
        //    	//printf(" coordD.z = %f\n", coordD.z);
        //        	printf("i.x = %d", i.x);
        //        	printf("  i.y = %d", i.y);
        //        	printf("  i.z = %d\n", i.z);
        //        	//printf("kk = %d\n", kk);
        //        	//printf("get_image_width(d_IImage) = %d\n", get_image_width(d_IImage));
        //        	//printf(" get_image_height(d_IImage) = %d\n", get_image_height(d_IImage));
        //}
        if (apu != 0.f) {
#ifdef MEANDISTANCEBP
            apu = (apu + d_meanV[kk] * area);
            //temp += (apu + d_meanV[kk] * area * convert_float(get_image_width(d_IImage) * get_image_height(d_IImage))) * kerroin;
#endif
            temp += apu * kerroin;
            if (no_norm == 0u)
                wSum += kerroin * area;
        }
        //pz += nZ;
    }
    d_forw[idx] = temp;
    if (no_norm == 0u)
        d_Summ[idx] = wSum;
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