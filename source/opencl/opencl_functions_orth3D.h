
/*******************************************************************************************************************************************
* Special functions for the 3D orthogonal distance-based ray tracer and volume of intersection ray tracer.
*
* Copyright (C) 2020-2026 Ville-Veikko Wettenhovi, Niilo Saarlemo
*
* This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************************************************************************/

#define THR 0.01f

#define SQRT_8LN2 2.35482004503094938f // 2*sqrt(2*ln(2)), i.e. FWHM/sigma
#define INV_8LN2 0.18033688011112042f // 1/(8*ln(2))

#ifndef PTYPE4
#ifdef SPECT
// In this function the ray is defined as v0+t*v1, where v0 is the source end of the ray and v0+v1 is the detector end of the ray. Thus 0<=t<=1. The parallel distance d is computed as the distance from the source to the projection of point p onto the ray. That is, d = v1*(p-v0) / (v1*v1) * |v1|
// return dot(v1, p - v0) * length(v1) / dot(v1, v1);
DEVICE FLOAT compute_element_parallel_3D(
	const FLOAT3 v0, // Ray begin point
	const FLOAT3 v1, // Ray end point - begin point
	const FLOAT3 p, // Voxel centre
	const FLOAT orth_ray_length_inv_signed // length(v1) / dot(v1, v1)
) {
    return dot(v1, p - v0) * orth_ray_length_inv_signed;
}

#define INV_SQRT_2PI 0.3989422804014327f // 1/sqrt(2*pi)
#define INV_2PI 0.15915494309189535f // 1/(2*pi) = (1/sqrt(2*pi))^2
#define INV_2SQRT2LN2 0.42466090014400953f // 1/(2*sqrt(2*ln2))
DEVICE float normPDF2(const float x, const float mu, const float sigma) {
    const float inv_sigma = RCP(sigma);
    const float inv_sigma2 = inv_sigma * inv_sigma;
    const float a = (x - mu) * inv_sigma;
    const float a2 = a * a;
    const float e = EXP(-0.5f * a2);
    return INV_2PI * inv_sigma2 * e;
}

#if defined(ORTH) && !defined(VOL)
DEVICE bool spectOrthRayBoxInterval(
    const FLOAT3 p0,
    const FLOAT3 dir,
    const FLOAT3 boxMin,
    const FLOAT3 boxMax,
    PTR_THR FLOAT *tmin,
    PTR_THR FLOAT *tmax
) {
    *tmin = FLOAT_ZERO;
    *tmax = 1.0e8f;
    const FLOAT epsVal = 1.0e-6f;

    if (FABS(dir.x) < epsVal) {
        if (p0.x < boxMin.x || p0.x > boxMax.x)
            return false;
    } else {
        const FLOAT t1 = (boxMin.x - p0.x) / dir.x;
        const FLOAT t2 = (boxMax.x - p0.x) / dir.x;
        *tmin = FMAX(*tmin, FMIN(t1, t2));
        *tmax = FMIN(*tmax, FMAX(t1, t2));
        if (*tmax < *tmin)
            return false;
    }

    if (FABS(dir.y) < epsVal) {
        if (p0.y < boxMin.y || p0.y > boxMax.y)
            return false;
    } else {
        const FLOAT t1 = (boxMin.y - p0.y) / dir.y;
        const FLOAT t2 = (boxMax.y - p0.y) / dir.y;
        *tmin = FMAX(*tmin, FMIN(t1, t2));
        *tmax = FMIN(*tmax, FMAX(t1, t2));
        if (*tmax < *tmin)
            return false;
    }

    if (FABS(dir.z) < epsVal) {
        if (p0.z < boxMin.z || p0.z > boxMax.z)
            return false;
    } else {
        const FLOAT t1 = (boxMin.z - p0.z) / dir.z;
        const FLOAT t2 = (boxMax.z - p0.z) / dir.z;
        *tmin = FMAX(*tmin, FMIN(t1, t2));
        *tmax = FMIN(*tmax, FMAX(t1, t2));
        if (*tmax < *tmin)
            return false;
    }

    return true;
}

DEVICE FLOAT spectOrthSupportRadius(
    const FLOAT3 s,
    const FLOAT3 diff,
    const FLOAT3 boxMin,
    const FLOAT3 boxMax,
    const FLOAT coneOfResponseStdCoeffA,
    const FLOAT coneOfResponseStdCoeffB,
    const FLOAT coneOfResponseStdCoeffC
) {
    const FLOAT rayLength = LENGTH(diff);
    if (rayLength <= 1.0e-6f)
        return FLOAT_ZERO;

    const FLOAT3 dir = diff * RCP(rayLength);
    FLOAT dParallelMax = FLOAT_ZERO;

    for (int ix = 0; ix < 2; ix++) {
        const FLOAT x = (ix == 0) ? boxMin.x : boxMax.x;
        for (int iy = 0; iy < 2; iy++) {
            const FLOAT y = (iy == 0) ? boxMin.y : boxMax.y;
            for (int iz = 0; iz < 2; iz++) {
                const FLOAT z = (iz == 0) ? boxMin.z : boxMax.z;
                const FLOAT3 corner = CMFLOAT3(x, y, z);
                dParallelMax = FMAX(dParallelMax, dot(dir, corner - s));
            }
        }
    }

    const FLOAT fwhm = FMAD(coneOfResponseStdCoeffA, dParallelMax, coneOfResponseStdCoeffB);
    return 3.5f * SQRT(fwhm * fwhm + coneOfResponseStdCoeffC * coneOfResponseStdCoeffC) * INV_2SQRT2LN2;
}

DEVICE bool spectOrthPrimaryIsX(const FLOAT3 diff) {
    return FABS(diff.x) >= FABS(diff.y);
}

DEVICE int spectOrthClampedIndex(const FLOAT coord, const FLOAT b0, const FLOAT d0, const uint n0) {
    const int maxIndex = CINT(n0) - 1;
    return MIN(MAX(CINT_rtz((coord - b0) / d0), 0), maxIndex);
}
#endif

#if defined(ODRT_TEXTURE)
#define SPECTPSF_RESOURCE IMAGE3D
#else
#define SPECTPSF_RESOURCE const CLGLOBAL float* CLRESTRICT
#endif

// Zero-extended read from a Fortran/ArrayFire ordered PSF volume. The first
// axis is contiguous: u + Nu * (v + Nv * depth).
#if !defined(ODRT_TEXTURE)
DEVICE float readSPECTPSF3DZero(const CLGLOBAL float* psf,
	const uint Nu, const uint Nv, const uint Nd,
	const int iu, const int iv, const int id) {
	if (psf == 0 || iu < 0 || iv < 0 || id < 0 || iu >= CINT(Nu) || iv >= CINT(Nv) || id >= CINT(Nd))
		return 0.f;
	const size_t ind = (size_t)iu + (size_t)Nu * ((size_t)iv + (size_t)Nv * (size_t)id);
	return psf[ind];
}
#endif

// Sample a SPECT PSF in the ray-local orthogonal plane and along-ray depth.
// u/v are physical signed offsets in mm, centered at the lateral PSF center
// ((N - 1) / 2); depth is the physical distance from the original SPECT ray
// start to the voxel. Ellipse clipping advances the integration ray start, so
// the per-ray offset is added back before sampling. du/dv/dd are the sample
// pitches in mm. Trilinear interpolation uses zero extension:
// each out-of-range corner contributes zero, including within one sample pitch
// of the table boundary. This gives a smooth falloff at finite PSF support.
DEVICE float sampleSPECTPSF3D(SPECTPSF_RESOURCE psf,
	const uint Nu, const uint Nv, const uint Nd,
	const float u_mm, const float v_mm, const float depth_mm,
	const float du_mm, const float dv_mm, const float dd_mm) {
	if (Nu == 0u || Nv == 0u || Nd == 0u ||
		!(du_mm > 0.f) || !(dv_mm > 0.f) || !(dd_mm > 0.f))
		return 0.f;

	const float fu = u_mm / du_mm + 0.5f * ((float)Nu - 1.f);
	const float fv = v_mm / dv_mm + 0.5f * ((float)Nv - 1.f);
	const float fd = depth_mm / dd_mm;
	// If the coordinate is farther than one sample beyond the table, all
	// interpolation corners are zero. These comparisons also reject NaNs.
	if (!(fu >= -1.f && fu < (float)Nu) ||
		!(fv >= -1.f && fv < (float)Nv) ||
		!(fd >= -1.f && fd < (float)Nd))
		return 0.f;

#if defined(ODRT_TEXTURE)
	// CUDA's unnormalized texture coordinates place texel centers at n+0.5.
	// Border addressing supplies zero-valued neighbors for edge interpolation.
	return tex3D<float>(psf, fu + 0.5f, fv + 0.5f, fd + 0.5f);
#else

	const float bu = FLOOR(fu), bv = FLOOR(fv), bd = FLOOR(fd);
	const int iu0 = CINT(bu), iv0 = CINT(bv), id0 = CINT(bd);
	const int iu1 = iu0 + 1, iv1 = iv0 + 1, id1 = id0 + 1;
	const float tu = fu - bu, tv = fv - bv, td = fd - bd;

	const float p000 = readSPECTPSF3DZero(psf, Nu, Nv, Nd, iu0, iv0, id0);
	const float p100 = readSPECTPSF3DZero(psf, Nu, Nv, Nd, iu1, iv0, id0);
	const float p010 = readSPECTPSF3DZero(psf, Nu, Nv, Nd, iu0, iv1, id0);
	const float p110 = readSPECTPSF3DZero(psf, Nu, Nv, Nd, iu1, iv1, id0);
	const float p001 = readSPECTPSF3DZero(psf, Nu, Nv, Nd, iu0, iv0, id1);
	const float p101 = readSPECTPSF3DZero(psf, Nu, Nv, Nd, iu1, iv0, id1);
	const float p011 = readSPECTPSF3DZero(psf, Nu, Nv, Nd, iu0, iv1, id1);
	const float p111 = readSPECTPSF3DZero(psf, Nu, Nv, Nd, iu1, iv1, id1);

	const float p00 = p000 + tu * (p100 - p000);
	const float p10 = p010 + tu * (p110 - p010);
	const float p01 = p001 + tu * (p101 - p001);
	const float p11 = p011 + tu * (p111 - p011);
	const float p0 = p00 + tv * (p10 - p00);
	const float p1 = p01 + tv * (p11 - p01);
	return p0 + td * (p1 - p0);
#endif
}

// Intersect an integer scan interval with the set of indices whose affine
// coordinate lies within one table axis' zero-extended interpolation support.
// Extra indices on both sides guard float rounding and keep this independent
// of the actual PSF values (which may contain holes).
DEVICE bool clipSPECTPSFInterval(const FLOAT coord0, const FLOAT coordStep,
	const FLOAT lower, const FLOAT upper, PTR_THR int* lo, PTR_THR int* hi) {
	if (*lo >= *hi)
		return false;
	if (FABS(coordStep) < 1.0e-8f) {
		if (coord0 < lower - 1.0e-4f || coord0 > upper + 1.0e-4f) {
			*hi = *lo;
			return false;
		}
		return true;
	}
	const FLOAT t0 = (lower - coord0) / coordStep;
	const FLOAT t1 = (upper - coord0) / coordStep;
	const FLOAT tlo = FMIN(t0, t1);
	const FLOAT thi = FMAX(t0, t1);
	const int clippedLo = MAX(*lo, CINT(CEIL(tlo)) - 1);
	const int clippedHi = MIN(*hi, CINT(FLOOR(thi)) + 2);
	*lo = clippedLo;
	*hi = MAX(clippedLo, clippedHi);
	return *lo < *hi;
}
#endif

// compute voxel index, orthogonal distance based or volume of intersection ray tracer
DEVICE LONG compute_ind_orth_3D(const uint tempi, const uint tempijk, const int tempk, const uint d_N, const uint Nyx) {
	LONG local_ind = CLONG_rtz(tempi * d_N + tempijk) + CLONG_rtz(tempk) * CLONG_rtz(Nyx);
	return local_ind;
}

#if !defined(SPECT) && defined(CRYSTXY)
// This function determines the number of voxels that are inside the "tube" in each slice
// This is done by determining the intersection points of the line and the cylinder
// These points are then converted into voxel indices
// Note that the "line" here is the line along the voxel axis investigated, not the actual ray
// Previously this was performed on the fly, where the distance from each voxel center to the "tube"
// was computed
DEVICE void tubeRangeXY(const FLOAT quadA, const FLOAT quadT, const FLOAT3 l, const FLOAT3 diff,
	const FLOAT d2, const FLOAT apuXY, const int minimiXY, const int maksimiXY, PTR_THR int* loXY, PTR_THR int* hiXY) {
	*loXY = minimiXY;
	*hiXY = maksimiXY;
	if (quadA <= FLOAT_ZERO) // Both diff.x and diff.z are zero
		return;
	// We use the distance to determine the intersection points and voxel indices
	// This requires solving a quadratic equation
	// These refer to the usual A, B, C letters in the quadratic formula
	const FLOAT quadB = -FLOAT_TWO * (diff.z * l.y + diff.x * l.z);
	const FLOAT quadC = l.x * l.x + l.y * l.y + l.z * l.z;
	// The discriminant
	// if this is negative, the line never intersects the "tube"
	// Zero is a special case of tangential line, but that should never occur in practice
	const FLOAT disc = quadB * quadB - 4.f * quadA * (quadC - quadT);
	if (disc <= FLOAT_ZERO) { // The "tube" does not reach this slice
		*hiXY = *loXY;
		return;
	}
	const FLOAT sq = SQRT(disc);
	const FLOAT inv2Ad = FLOAT_HALF / (quadA * d2);
	// Here the distance is converted into the voxel indices
#ifdef USEMAD
	const FLOAT loF = FMAX(FMAD((-quadB - sq), inv2Ad, apuXY), CFLOAT(minimiXY));
	const FLOAT hiF = FMIN(FMAD((-quadB + sq), inv2Ad, apuXY), CFLOAT(maksimiXY));
#else
	const FLOAT loF = FMAX((-quadB - sq) * inv2Ad + apuXY, CFLOAT(minimiXY));
	const FLOAT hiF = FMIN((-quadB + sq) * inv2Ad + apuXY, CFLOAT(maksimiXY));
#endif
	// Clamp the indices
	*loXY = MAX(CINT(CEIL(loF)), minimiXY);
	*hiXY = MIN(CINT(FLOOR(hiF)) + 1, maksimiXY);
	if (*hiXY < *loXY)
		*hiXY = *loXY;
}
#endif

// This function computes either the forward or backward projection for the current voxel
// The normalized orthogonal distance or the volume of the (spherical) voxel is computed before the forward or backward projections
DEVICE bool orthogonalHelper3D(const int tempi, const int uu, const uint d_N2, const uint d_N3, const uint d_Nxy, const int zz, 
    const FLOAT3 s,
    const FLOAT3 l,
	const FLOAT3 diff, 
    const FLOAT orth_ray_length, const FLOAT orth_ray_length_inv, const FLOAT3 center, const FLOAT bmin, const FLOAT bmax, const FLOAT Vmax, CONSTANT float* V, const bool XY, PTR_THR float *ax, const FLOAT temp, 
#if defined(FP)
	IMTYPE d_OSEM
#else
	CLGLOBAL CAST* d_Summ, CLGLOBAL CAST* d_output, const bool no_norm
#endif
#ifdef TOF
	, const float TOFSum, float* TOFWeights
#ifdef LISTMODE
	, const int TOFid
#endif
#endif
#if defined(SPECT) && defined(ATN)
	, FLOAT jelppi
#endif
#if defined(MASKBP) && defined(BP)
    , const int ii, MASKBPTYPE maskBP, const uint3 d_N
#endif
#ifdef SPECT
    , const FLOAT coneOfResponseStdCoeffA, const FLOAT coneOfResponseStdCoeffB, const FLOAT c2, const FLOAT orth_ray_length_inv_signed
	, const FLOAT rayDepthOffset
	, SPECTPSF_RESOURCE gFilter, const uint Nu, const uint Nv, const uint Nd
	, const FLOAT du, const FLOAT dv, const FLOAT dd, const bool gFilterCustom
	, const FLOAT3 ray, const FLOAT3 e1, const FLOAT3 e2
#endif
) {
#ifdef SPECT // Check for voxel behind detector
    float d_parallel = compute_element_parallel_3D(s, diff, center, orth_ray_length_inv_signed);
    if (d_parallel < 0) { 
        return false;
    }
#endif
	float local_ele;
#ifdef SPECT ////////////////////////// SPECT ////////////////////
	if (gFilterCustom) {
		const FLOAT3 q = (center - s) - ray * d_parallel;
		const float u = dot(q, e1);
		const float v = dot(q, e2);
		local_ele = sampleSPECTPSF3D(gFilter, Nu, Nv, Nd, u, v, d_parallel + rayDepthOffset, du, dv, dd);
	} else {
		const FLOAT x0 = FMAD(-FLOAT_ONE, s.y, center.y);
		const FLOAT y1 = FMAD(diff.z, x0, -l.y);
		const FLOAT z1 = FMAD(-diff.x, x0, l.z);
		const FLOAT norm2 = FMAD(l.x, l.x, FMAD(y1, y1, z1 * z1));
		const float t = FMAD(coneOfResponseStdCoeffA, d_parallel, coneOfResponseStdCoeffB);   // A*d_parallel + B
		const float var8ln2 = FMAD(t, t, c2); // (2*sqrt(2*ln2)*CORstd)^2
		// norm2 >= 3.5^2 * CORstd^2 * orth_ray_length^2;
		if (norm2 >= 12.25f * INV_8LN2 * var8ln2 * orth_ray_length * orth_ray_length) {
			return true;
		}
		const float invSTD = RSQRT(var8ln2) * SQRT_8LN2; // 1/CORstd
		const float a = SQRT(norm2) * orth_ray_length_inv * invSTD;
		local_ele = INV_2PI * invSTD * invSTD * EXP(-0.5f * a * a);
	}
#elif defined(VOL) //////////// VOL /////////////
	const FLOAT x0 = FMAD(-FLOAT_ONE, s.y, center.y);
	const FLOAT y1 = FMAD(diff.z, x0, -l.y);
	const FLOAT z1 = FMAD(-diff.x, x0, l.z);
	const FLOAT norm2 = FMAD(l.x, l.x, FMAD(y1, y1, z1 * z1));
	if (norm2 > bmax * bmax * orth_ray_length * orth_ray_length) {
		return true;
	}
	local_ele = SQRT(norm2) * orth_ray_length_inv;
	if (local_ele < bmin)
		local_ele = Vmax;
	else
		local_ele = V[CUINT_rte((FMIN(local_ele, bmax) - bmin) * CC)];
#else //////////// ORTH /////////////
	const FLOAT x0 = FMAD(-FLOAT_ONE, s.y, center.y);
	const FLOAT y1 = FMAD(diff.z, x0, -l.y);
	const FLOAT z1 = FMAD(-diff.x, x0, l.z);
	const FLOAT norm2 = FMAD(l.x, l.x, FMAD(y1, y1, z1 * z1));
	if (norm2 >= (FLOAT_ONE - THR) * (FLOAT_ONE - THR) * orth_ray_length * orth_ray_length) {
		return true;
	}
	local_ele = FLOAT_ONE - SQRT(norm2) * orth_ray_length_inv;
#endif //////////// END SPECT/VOL/ORTH /////////////
#if (defined(FP) && defined(USEIMAGES)) || (defined(MASKBP) && defined(BP)) ///////////////////// 2D/3D indices /////////////////////
	int3 ind;
    if (XY)
        ind = CMINT3(tempi, uu, zz);
    else
        ind = CMINT3(uu, tempi, zz);
#endif
#if defined(BP) || (defined(FP) && !defined(USEIMAGES))
	LONG local_ind = 0;
	local_ind = compute_ind_orth_3D(CUINT(tempi), uu * d_N3, (zz), d_N2, d_Nxy);
#endif
#if (defined(FP) && !defined(USEIMAGES))
	const LONG ind = CLONG_rtz(local_ind);
#endif ///////////////////// END 2D/3D indices /////////////////////
#if defined(BP) && defined(MASKBP) ///////////////////// APPLY BP MASK /////////////////////////
    if (readMaskBP(maskBP, ind, d_N) == 0) {
        return false;
    }
#endif ///////////////////// END BP MASK /////////////////////////
#if defined(SPECT) && defined(ATN)
	local_ele *= jelppi;
#endif
#if defined(FP) //////////////// FORWARD PROJECTION ////////////////
	denominator(ax, ind, local_ele, d_OSEM
#ifdef TOF //////////////// TOF ////////////////
				, TOFSum, TOFWeights
#ifdef LISTMODE
				, TOFid
#endif
#endif //////////////// END TOF ////////////////
	);
#endif //////////////// END FORWARD PROJECTION ////////////////
#if defined(BP) //////////////// BACKWARD PROJECTION ////////////////
	rhs(local_ele * temp, ax, local_ind, d_output, no_norm, d_Summ
#ifdef TOF
			, TOFSum, TOFWeights
#ifdef LISTMODE
			, TOFid
#endif
#endif
    );
#endif //////////////// END BACKWARD PROJECTION ////////////////
	return false;
}
 
// Both the orthogonal and volume of intersection ray tracers loop through all the neighboring voxels of the current voxel
// Both also proceed through each X or Y slice, depending on the incident direction
// This function simply loops through each X or Y voxel and Z voxels in the current slice
// Forward or backward projection is computed in the helper function
DEVICE int orthDistance3D(const int tempi, 
    const FLOAT3 diff, // detectors.[x/y/z]d - detectors.[x/y/z]s
    FLOAT3 center, // center of the current voxel
    const FLOAT3 s, // detectors.[x/y/z]s
    const FLOAT b2, const FLOAT d2, const FLOAT bz, const FLOAT dz, const FLOAT temp, int temp2, const int tempk, 
	const uint d_Nxy, const FLOAT orth_ray_length,
	const uint d_N1, const uint d_N2, const uint d_N3, const uint d_Nz, const FLOAT bmin, 
	const FLOAT bmax, const FLOAT Vmax, CONSTANT float* V, const bool XY, PTR_THR float *ax, const bool preStep, PTR_THR int *k, const int ku, 
#if defined(FP)
	IMTYPE d_OSEM
#else
	const bool no_norm, CLGLOBAL CAST* Summ, CLGLOBAL CAST* d_rhs_OSEM 
#endif
#ifdef TOF
	, const float TOFSum, float* TOFWeights
#ifdef LISTMODE
	, const int TOFid
#endif
#endif
#if defined(SPECT) && defined(ATN)
	, FLOAT jelppi
#endif
#if defined(MASKBP) && defined(BP)
	, const int ii, MASKBPTYPE maskBP, const uint3 d_N
#endif
#ifdef SPECT
    , const FLOAT coneOfResponseStdCoeffA, const FLOAT coneOfResponseStdCoeffB, const FLOAT coneOfResponseStdCoeffC, const FLOAT orth_ray_length_inv_signed
	, const FLOAT rayDepthOffset
	, SPECTPSF_RESOURCE gFilter, const uint Nu, const uint Nv, const uint Nd
	, const FLOAT du, const FLOAT dv, const FLOAT dd, const bool gFilterCustom
#endif
) {
	int uu = 0;
	bool breikki = false;
#ifdef SPECT
	const FLOAT c2 = coneOfResponseStdCoeffC * coneOfResponseStdCoeffC;
#endif
    FLOAT3 l; // Precomputed cross product elements
	const FLOAT v0 = center.x - s.x;
    l.z = diff.y * v0;
	const FLOAT apu1 = diff.z * v0;
	const int maksimiZ = CINT(d_Nz);
	const int minimiZ = 0;
	const int maksimiXY = CINT(d_N1);
	const int minimiXY = 0;
	int uu1 = 0, uu2 = 0;
	// Computed once here rather than dividing separately for every voxel in the helper
	const FLOAT orth_ray_length_inv = FLOAT_ONE / orth_ray_length;
#ifdef SPECT
	FLOAT3 ray = MFLOAT3(FLOAT_ZERO, FLOAT_ZERO, FLOAT_ZERO);
	FLOAT3 e1 = MFLOAT3(FLOAT_ZERO, FLOAT_ZERO, FLOAT_ZERO);
	FLOAT3 e2 = MFLOAT3(FLOAT_ZERO, FLOAT_ZERO, FLOAT_ZERO);
	if (gFilterCustom) {
		ray = diff * orth_ray_length_inv;
		FLOAT3 ref = MFLOAT3(FLOAT_ZERO, FLOAT_ZERO, FLOAT_ONE);
		if (FABS(ray.z) > 0.999f)
			ref = MFLOAT3(FLOAT_ONE, FLOAT_ZERO, FLOAT_ZERO);
		e1 = NORMALIZE(CROSS(ref, ray));
		e2 = CROSS(ray, e1);
	}
#endif
	int loXY = minimiXY, hiXY = maksimiXY;
#if !defined(SPECT) && defined(CRYSTXY)
	const FLOAT quadA = diff.x * diff.x + diff.z * diff.z;
#if defined(VOL)
	const FLOAT quadT = bmax * bmax * orth_ray_length * orth_ray_length;
#else
	const FLOAT quadT = (FLOAT_ONE - THR) * (FLOAT_ONE - THR) * orth_ray_length * orth_ray_length;
#endif
	// The distance is solved for x0 = center.y - s.y, where center.y = b2 + uu * d2 + d2 / 2
	// that is, the center of the voxel
	// Independent of z, so it is computed once here rather than in every slice
	const FLOAT apuXY = (s.y - b2) / d2 - FLOAT_HALF;
#endif
	int zz = tempk;
#ifdef CRYSTZ
	for (zz = MAX(tempk, minimiZ); zz < maksimiZ; zz++) {
#endif
		center.z = bz + CFLOAT(zz) * dz + dz / FLOAT_TWO;
		int scanLoXY = loXY, scanHiXY = hiXY;
#if defined(SPECT) && defined(CRYSTXY)
		if (gFilterCustom) {
			// Each ray-local coordinate is affine in the transaxial scan index.
			// Restrict this slice to the table's zero-extended support; this is
			// independent of PSF values, including sparse/non-monotone tables.
			FLOAT3 origin = center;
			origin.y = b2 + d2 / FLOAT_TWO;
			const FLOAT3 offset = origin - s;
			const FLOAT parallel = dot(offset, ray);
			const FLOAT depth = parallel + rayDepthOffset;
			const FLOAT3 q = offset - ray * parallel;
			const FLOAT uLimit = FLOAT_HALF * (CFLOAT(Nu) + FLOAT_ONE) * du;
			const FLOAT vLimit = FLOAT_HALF * (CFLOAT(Nv) + FLOAT_ONE) * dv;
			const FLOAT depthLimit = CFLOAT(Nd) * dd;
			bool hasSupport = clipSPECTPSFInterval(dot(q, e1), d2 * e1.y,
				-uLimit, uLimit, &scanLoXY, &scanHiXY);
			if (hasSupport)
				hasSupport = clipSPECTPSFInterval(dot(q, e2), d2 * e2.y,
					-vLimit, vLimit, &scanLoXY, &scanHiXY);
			if (hasSupport)
				hasSupport = clipSPECTPSFInterval(depth, d2 * ray.y,
					FLOAT_ZERO, depthLimit, &scanLoXY, &scanHiXY);
			if (!hasSupport)
				scanHiXY = scanLoXY;
		}
#endif
		const FLOAT z0 = center.z - s.z;
		l.x = diff.x * z0 - apu1;
		l.y = diff.y * z0;
		// Compute the maximum ranges
#if !defined(SPECT) && defined(CRYSTXY)
		tubeRangeXY(quadA, quadT, l, diff, d2, apuXY, minimiXY, maksimiXY, &loXY, &hiXY);
#endif
#ifdef CRYSTXY
#if !defined(SPECT)
		// tubeRangeXY gives the exact range, so one ascending pass covers the whole row. The
		// bidirectional scan below is only needed when the tube end has to be found on the fly
		for (uu1 = scanLoXY; uu1 < scanHiXY; uu1++) {
#else
		for (uu1 = MAX(temp2, scanLoXY); uu1 < scanHiXY; uu1++) {
#endif
#else
		uu1 = temp2;
#endif
			center.y = b2 + CFLOAT(uu1) * d2 + d2 / FLOAT_TWO;
			breikki = orthogonalHelper3D(tempi, uu1, d_N2, d_N3, d_Nxy, zz, s, l, diff, orth_ray_length, orth_ray_length_inv, center, bmin, bmax, Vmax, V, XY, ax, temp, 
#if defined(FP)
				d_OSEM
#else
				Summ, d_rhs_OSEM, no_norm
#endif
#ifdef TOF
				, TOFSum, TOFWeights
#ifdef LISTMODE
				, TOFid
#endif
#endif
#if defined(SPECT) && defined(ATN)
				, jelppi
#endif
#if defined(MASKBP) && defined(BP)
				, ii, maskBP, d_N
#endif
#ifdef SPECT
				, coneOfResponseStdCoeffA, coneOfResponseStdCoeffB, c2, orth_ray_length_inv_signed
				, rayDepthOffset, gFilter, Nu, Nv, Nd, du, dv, dd, gFilterCustom
				, ray, e1, e2
#endif
			);
#ifdef CRYSTXY
#if !defined(SPECT)
			// tubeRangeXY already bounds this loop exactly, so a voxel rejected here can only be the
			// float-rounding boundary case. Breaking would drop the remainder of the row, so it is skipped
			if (!breikki)
				uu++;
#else
			if (breikki) {
				break;
			}
			uu++;
#endif
		}
#else
		// The caller stops scanning slices once the return value is zero, so the single voxel of the
		// no transaxial spread case has to be counted as well
		if (!breikki)
			uu++;
#endif
#if defined(CRYSTXY) && defined(SPECT)
		// Only SPECT reaches this: without the analytic range the scan has to start at the ray and
		// walk outward in both directions, stopping at the first voxel outside the response
		for (uu2 = MIN(temp2, scanHiXY) - 1; uu2 >= scanLoXY; uu2--) {
			center.y = b2 + CFLOAT(uu2) * d2 + d2 / FLOAT_TWO;
			breikki = orthogonalHelper3D(tempi, uu2, d_N2, d_N3, d_Nxy, zz, s, l, diff, orth_ray_length, orth_ray_length_inv, center, bmin, bmax, Vmax, V, XY, ax, temp, 
#if defined(FP)
				d_OSEM
#else
				Summ, d_rhs_OSEM, no_norm
#endif
#ifdef TOF
				, TOFSum, TOFWeights
#ifdef LISTMODE
				, TOFid
#endif
#endif
#if defined(SPECT) && defined(ATN)
				, jelppi
#endif
#if defined(MASKBP) && defined(BP)
				, ii, maskBP, d_N
#endif
#ifdef SPECT
				, coneOfResponseStdCoeffA, coneOfResponseStdCoeffB, c2, orth_ray_length_inv_signed
				, rayDepthOffset, gFilter, Nu, Nv, Nd, du, dv, dd, gFilterCustom
				, ray, e1, e2
#endif
			);
			if (breikki) {
				break;
			}
			uu++;
		}
#elif !defined(CRYSTXY)
	uu2 = temp2 - 1;
#endif
#ifdef CRYSTZ
#if !defined(SPECT) && defined(CRYSTXY)
	// The transaxial range is exact, so an empty range means the tube has ended in the axial direction too
	if (loXY >= hiXY)
		break;
#else
#if defined(SPECT) && defined(CRYSTXY)
		if (!gFilterCustom && uu1 == temp2 && uu2 == temp2 - 1 && breikki)
			break;
#else
		if (uu1 == temp2 && uu2 == temp2 - 1 && breikki)
			break;
#endif
#endif
	}
	*k = zz - 1;
	for (zz = MIN(tempk, maksimiZ) - 1; zz >= minimiZ; zz--) {
		center.z = bz + CFLOAT(zz) * dz + dz / FLOAT_TWO;
		int scanLoXY = loXY, scanHiXY = hiXY;
#if defined(SPECT) && defined(CRYSTXY)
		if (gFilterCustom) {
			FLOAT3 origin = center;
			origin.y = b2 + d2 / FLOAT_TWO;
			const FLOAT3 offset = origin - s;
			const FLOAT parallel = dot(offset, ray);
			const FLOAT depth = parallel + rayDepthOffset;
			const FLOAT3 q = offset - ray * parallel;
			const FLOAT uLimit = FLOAT_HALF * (CFLOAT(Nu) + FLOAT_ONE) * du;
			const FLOAT vLimit = FLOAT_HALF * (CFLOAT(Nv) + FLOAT_ONE) * dv;
			const FLOAT depthLimit = CFLOAT(Nd) * dd;
			bool hasSupport = clipSPECTPSFInterval(dot(q, e1), d2 * e1.y,
				-uLimit, uLimit, &scanLoXY, &scanHiXY);
			if (hasSupport)
				hasSupport = clipSPECTPSFInterval(dot(q, e2), d2 * e2.y,
					-vLimit, vLimit, &scanLoXY, &scanHiXY);
			if (hasSupport)
				hasSupport = clipSPECTPSFInterval(depth, d2 * ray.y,
					FLOAT_ZERO, depthLimit, &scanLoXY, &scanHiXY);
			if (!hasSupport)
				scanHiXY = scanLoXY;
		}
#endif
		const FLOAT z0 = center.z - s.z;
		l.x = diff.x * z0 - apu1;
		l.y = diff.y * z0;
#if !defined(SPECT) && defined(CRYSTXY)
		tubeRangeXY(quadA, quadT, l, diff, d2, apuXY, minimiXY, maksimiXY, &loXY, &hiXY);
#endif
#ifdef CRYSTXY
#if !defined(SPECT)
		// tubeRangeXY gives the exact range, so one ascending pass covers the whole row. The
		// bidirectional scan below is only needed when the tube end has to be found on the fly
		for (uu1 = scanLoXY; uu1 < scanHiXY; uu1++) {
#else
		for (uu1 = MAX(temp2, scanLoXY); uu1 < scanHiXY; uu1++) {
#endif
#else
		uu1 = temp2;
#endif
			center.y = b2 + CFLOAT(uu1) * d2 + d2 / FLOAT_TWO;
			breikki = orthogonalHelper3D(tempi, uu1, d_N2, d_N3, d_Nxy, zz, s, l, diff, orth_ray_length, orth_ray_length_inv, center, bmin, bmax, Vmax, V, XY, ax, temp, 
#if defined(FP)
				d_OSEM
#else
				Summ, d_rhs_OSEM, no_norm
#endif
#ifdef TOF
				, TOFSum, TOFWeights
#ifdef LISTMODE
				, TOFid
#endif
#endif
#if defined(SPECT) && defined(ATN)
				, jelppi
#endif
#if defined(MASKBP) && defined(BP)
				, ii, maskBP, d_N
#endif
#ifdef SPECT
				, coneOfResponseStdCoeffA, coneOfResponseStdCoeffB, c2, orth_ray_length_inv_signed
				, rayDepthOffset, gFilter, Nu, Nv, Nd, du, dv, dd, gFilterCustom
				, ray, e1, e2
#endif
			);
#ifdef CRYSTXY
#if !defined(SPECT)
			// tubeRangeXY already bounds this loop exactly, so a voxel rejected here can only be the
			// float-rounding boundary case. Breaking would drop the remainder of the row, so it is skipped
			if (!breikki)
				uu++;
#else
			if (breikki) {
				break;
			}
			uu++;
#endif
		}
#else
		// The caller stops scanning slices once the return value is zero, so the single voxel of the
		// no transaxial spread case has to be counted as well
		if (!breikki)
			uu++;
#endif
#if defined(CRYSTXY) && defined(SPECT)
		// Only SPECT reaches this: without the analytic range the scan has to start at the ray and
		// walk outward in both directions, stopping at the first voxel outside the response
		for (uu2 = MIN(temp2, scanHiXY) - 1; uu2 >= scanLoXY; uu2--) {
			center.y = b2 + CFLOAT(uu2) * d2 + d2 / FLOAT_TWO;
			breikki = orthogonalHelper3D(tempi, uu2, d_N2, d_N3, d_Nxy, zz, s, l, diff, orth_ray_length, orth_ray_length_inv, center, bmin, bmax, Vmax, V, XY, ax, temp, 
#if defined(FP)
				d_OSEM
#else
				Summ, d_rhs_OSEM, no_norm
#endif
#ifdef TOF
				, TOFSum, TOFWeights
#ifdef LISTMODE
				, TOFid
#endif
#endif
#if defined(SPECT) && defined(ATN)
				, jelppi
#endif
#if defined(MASKBP) && defined(BP)
				, ii, maskBP, d_N
#endif
#ifdef SPECT
				, coneOfResponseStdCoeffA, coneOfResponseStdCoeffB, c2, orth_ray_length_inv_signed
				, rayDepthOffset, gFilter, Nu, Nv, Nd, du, dv, dd, gFilterCustom
				, ray, e1, e2
#endif
			);
			if (breikki) {
				break;
			}
			uu++;
		}
#elif !defined(CRYSTXY)
	uu2 = temp2 - 1;
#endif
#if !defined(SPECT) && defined(CRYSTXY)
		if (loXY >= hiXY)
			break;
#else
	#if defined(SPECT) && defined(CRYSTXY)
		if (!gFilterCustom && uu1 == temp2 && uu2 == temp2 - 1 && breikki)
			break;
	#else
		if (uu1 == temp2 && uu2 == temp2 - 1 && breikki)
			break;
	#endif
#endif
	}
	if (preStep) {
		if (ku < 0)
			*k = MAX(*k, zz) - 1;
		else
			*k = MIN(*k, zz) + 1;
	}
#endif
	return uu;
}

#endif
