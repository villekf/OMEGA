
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

#ifdef SPECT // Threshold for ODRT collimator response / tube width cutoff
#define THR normPDF2(3.5*CORstd, 0.f, CORstd)
#else
#define THR 0.01f
#endif

#define SQRT_8LN2 2.35482004503094938f // 2*sqrt(2*ln(2)), i.e. FWHM/sigma
#define INV_8LN2 0.18033688011112042f // 1/(8*ln(2))

#ifndef PTYPE4
// Compute the orthogonal distance from the ray to the current voxel (center)
// For orthogonal distance-based ray tracer, the distance is normalized
// See for example https://math.stackexchange.com/questions/2353288/point-to-line-distance-in-3d-using-cross-product/2355960#2355960
// In this function the ray is defined as v0+t*v1, where v0 is the source end of the ray and v0+v1 is the detector end of the ray. Thus 0<=t<=1. The orthogonal distance d from the voxel center p to the ray is equal to d = |v1 x (v0-p)| / |v1|. The term |v1| is precomputed, and in PET reconstruction includes the crystal size.
// Precomputed:
//l.x = diff.x * (center.z - s.z) - diff.z  * (center.x - s.x);
//l.y = diff.y * (center.z - s.z);
//l.z = diff.y * (center.x - s.x);
DEVICE FLOAT compute_element_orth_3D(
    const FLOAT3 s, // Source position
    const FLOAT3 l, // Precomputed cross product elements
    const FLOAT3 diff, // Spans the ray
    const FLOAT3 center, // Center of the voxel
    const FLOAT orth_ray_length // Precomputed length |diff|, in PET scaled with tube width
) {
    const FLOAT x0 = FMAD(-FLOAT_ONE, s.y, center.y);
	const FLOAT y1 = FMAD(diff.z, x0, -l.y);
	const FLOAT z1 = FMAD(-diff.x, x0, l.z);
	const FLOAT norm1 = LENGTH(CMFLOAT3(l.x, y1, z1));
    return
#if !defined(VOL) && !defined(SPECT)
    FLOAT_ONE - 
#endif
    norm1 / orth_ray_length;
}

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
#endif

// compute voxel index, orthogonal distance based or volume of intersection ray tracer
DEVICE LONG compute_ind_orth_3D(const uint tempi, const uint tempijk, const int tempk, const uint d_N, const uint Nyx) {
	LONG local_ind = CLONG_rtz(tempi * d_N + tempijk) + CLONG_rtz(tempk) * CLONG_rtz(Nyx);
	return local_ind;
}

// This function computes either the forward or backward projection for the current voxel
// The normalized orthogonal distance or the volume of the (spherical) voxel is computed before the forward or backward projections
DEVICE bool orthogonalHelper3D(const int tempi, const int uu, const uint d_N2, const uint d_N3, const uint d_Nxy, const int zz, 
    const FLOAT3 s,
    const FLOAT3 l,
	const FLOAT3 diff, 
    const FLOAT orth_ray_length, const FLOAT3 center, const FLOAT bmin, const FLOAT bmax, const FLOAT Vmax, CONSTANT float* V, const bool XY, PTR_THR float *ax, const FLOAT temp, 
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
#endif
) {
#ifdef SPECT // Check for voxel behind detector
    float d_parallel = compute_element_parallel_3D(s, diff, center, orth_ray_length_inv_signed);
    if (d_parallel < 0) { 
        return false;
    }
#endif
    const FLOAT x0 = FMAD(-FLOAT_ONE, s.y, center.y);
	const FLOAT y1 = FMAD(diff.z, x0, -l.y);
	const FLOAT z1 = FMAD(-diff.x, x0, l.z);
	const FLOAT norm2 = FMAD(l.x, l.x, FMAD(y1, y1, z1 * z1));
	float local_ele;
#ifdef SPECT ////////////////////////// SPECT ////////////////////
	const float t = FMAD(coneOfResponseStdCoeffA, d_parallel, coneOfResponseStdCoeffB);   // A*d_parallel + B
	const float var8ln2 = FMAD(t, t, c2); // (2*sqrt(2*ln2)*CORstd)^2
	// norm2 >= 3.5^2 * CORstd^2 * orth_ray_length^2;
	if (norm2 >= 12.25f * INV_8LN2 * var8ln2 * orth_ray_length * orth_ray_length) {
		return true;
	}
	const float invSTD = RSQRT(var8ln2) * SQRT_8LN2; // 1/CORstd
	const float a = SQRT(norm2) / orth_ray_length * invSTD;
	local_ele = INV_2PI * invSTD * invSTD * EXP(-0.5f * a * a);
#elif defined(VOL) //////////// VOL /////////////
	if (norm2 > bmax * bmax * orth_ray_length * orth_ray_length) {
		return true;
	}
	local_ele = SQRT(norm2) / orth_ray_length;
	if (local_ele < bmin)
		local_ele = Vmax;
	else
		local_ele = V[CUINT_rte((local_ele - bmin) * CC)];
#else //////////// ORTH /////////////
	if (norm2 >= (FLOAT_ONE - THR) * (FLOAT_ONE - THR) * orth_ray_length * orth_ray_length) {
		return true;
	}
	local_ele = FLOAT_ONE - SQRT(norm2) / orth_ray_length;
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
    if ((ii == 0) && (readMaskBP(maskBP, ind, d_N) == 0)) {
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
	const FLOAT centerY0 = b2 + CFLOAT(temp2) * d2 + d2 / FLOAT_TWO;
	const FLOAT centerZ0 = bz + CFLOAT(tempk) * dz + dz / FLOAT_TWO;
	int zz = tempk;
	center.z = centerZ0;
#ifdef CRYSTZ
	for (zz = MAX(tempk, minimiZ); zz < maksimiZ; zz++) {
#endif
		center.z = bz + CFLOAT(zz) * dz + dz / FLOAT_TWO;
		const FLOAT z0 = center.z - s.z;
		l.x = diff.x * z0 - apu1;
		l.y = diff.y * z0;
#ifdef CRYSTXY
		center.y = centerY0;
		for (uu1 = MAX(temp2, minimiXY); uu1 < maksimiXY; uu1++) {
#else
		uu1 = temp2;
		center.y = centerY0;
#endif
			center.y = b2 + CFLOAT(uu1) * d2 + d2 / FLOAT_TWO;
			breikki = orthogonalHelper3D(tempi, uu1, d_N2, d_N3, d_Nxy, zz, s, l, diff, orth_ray_length, center, bmin, bmax, Vmax, V, XY, ax, temp, 
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
#endif
			);
#ifdef CRYSTXY
			if (breikki) {
				break;
			}
			uu++;
			center.y += d2;
		}
#endif
#ifdef CRYSTXY
		center.y = centerY0 - d2;
		for (uu2 = MIN(temp2, maksimiXY) - 1; uu2 >= minimiXY; uu2--) {
			breikki = orthogonalHelper3D(tempi, uu2, d_N2, d_N3, d_Nxy, zz, s, l, diff, orth_ray_length, center, bmin, bmax, Vmax, V, XY, ax, temp, 
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
#endif
			);
			if (breikki) {
				break;
			}
			uu++;
			center.y -= d2;
		}
#else
	uu2 = temp2 - 1;
#endif
#ifdef CRYSTZ
	if (uu1 == temp2 && uu2 == temp2 - 1 && breikki)
		break;
	center.z += dz;
	}
	*k = zz - 1;
	center.z = centerZ0 - dz;
	for (zz = MIN(tempk, maksimiZ) - 1; zz >= minimiZ; zz--) {
		const FLOAT z0 = center.z - s.z;
		l.x = diff.x * z0 - apu1;
		l.y = diff.y * z0;
#ifdef CRYSTXY
		center.y = centerY0;
		for (uu1 = MAX(temp2, minimiXY); uu1 < maksimiXY; uu1++) {
#else
		uu1 = temp2;
		center.y = centerY0;
#endif
			breikki = orthogonalHelper3D(tempi, uu1, d_N2, d_N3, d_Nxy, zz, s, l, diff, orth_ray_length, center, bmin, bmax, Vmax, V, XY, ax, temp, 
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
#endif
			);
#ifdef CRYSTXY
			if (breikki) {
				break;
			}
			uu++;
			center.y += d2;
		}
#endif
#ifdef CRYSTXY
		center.y = centerY0 - d2;
		for (uu2 = MIN(temp2, maksimiXY) - 1; uu2 >= minimiXY; uu2--) {
			breikki = orthogonalHelper3D(tempi, uu2, d_N2, d_N3, d_Nxy, zz, s, l, diff, orth_ray_length, center, bmin, bmax, Vmax, V, XY, ax, temp, 
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
#endif
			);
			if (breikki) {
				break;
			}
			uu++;
			center.y -= d2;
		}
#else
	uu2 = temp2 - 1;
#endif
		if (uu1 == temp2 && uu2 == temp2 - 1 && breikki)
			break;
		center.z -= dz;
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
