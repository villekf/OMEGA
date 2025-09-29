
/*******************************************************************************************************************************************
* Special functions for the 3D orthogonal distance-based ray tracer and volume of intersection ray tracer.
*
* Copyright (C) 2020-2024 Ville-Veikko Wettenhovi, Niilo Saarlemo
*
* This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************************************************************************/

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

#ifndef _2PI // 1/sqrt(2*pi)
#ifdef HALF
#define _2PI 0.3989423h
#else
#define _2PI 0.3989423f 
#endif
#endif
// Evaluate Gaussian distribution at x. Normalized by the 2D gaussian constant.
DEVICE FLOAT normPDF2(const FLOAT x, const FLOAT mu, const FLOAT sigma) {
    const FLOAT a = FMAD(-FLOAT_ONE, DIVIDE(mu, sigma), DIVIDE(x, sigma));
    return POWR(DIVIDE(_2PI, sigma), FLOAT_TWO) * EXP(-FLOAT_HALF * a * a);
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
	, const FLOAT element, const FLOAT sigma_x, float* D, const FLOAT DD, const FLOAT TOFSum, float* TOFWeights
#ifdef LISTMODE
	, const int TOFid
#endif
#endif
#if defined(SPECT) && defined(ATN)
	, FLOAT jelppi
#endif
#if defined(MASKBP) && defined(BP)
    , const int ii, MASKBPTYPE maskBP
#endif
#ifdef SPECT
    , const FLOAT coneOfResponseStdCoeffA, const FLOAT coneOfResponseStdCoeffB, const FLOAT coneOfResponseStdCoeffC, const FLOAT orth_ray_length_inv_signed
#endif
#ifdef N_RAYS
    , int lor
#endif
) {
#if (defined(FP) || (defined(MASKBP) && defined(BP))) && defined(USEIMAGES) ///////////////////// 2D/3D indices /////////////////////
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
#endif
#if (defined(BP) && defined(MASKBP)) && !defined(USEIMAGES)
    const LONG ind = CLONG_rtz(local_ind
#if !defined(MASKBP3D)
    - zz * d_Nxy
#endif
    );
#endif ///////////////////// END 2D/3D indices /////////////////////
	FLOAT local_ele = compute_element_orth_3D(s, l, diff, center, orth_ray_length);
#ifdef SPECT ////////////////////////// SPECT ////////////////////
    FLOAT d_parallel = compute_element_parallel_3D(s, diff, center, orth_ray_length_inv_signed);
    if (d_parallel < 0) { // Voxel behind detector
        return false;
    }
    FLOAT CORstd = DIVIDE(SQRT(FMAD(FLOAT_ONE, POWR(FMAD(coneOfResponseStdCoeffA, d_parallel, coneOfResponseStdCoeffB), FLOAT_TWO), POWR(coneOfResponseStdCoeffC, FLOAT_TWO))), (FLOAT_TWO*SQRT(FLOAT_TWO*LOG(FLOAT_TWO))));
    local_ele = normPDF2(local_ele, FLOAT_ZERO, CORstd);
#endif //////////// END SPECT /////////////
#ifdef VOL
	if (local_ele > bmax) {
		return true;
	}
	if (local_ele < bmin)
		local_ele = Vmax;
	else
		local_ele = V[CUINT_rte((local_ele - bmin) * CC)];
#else
#ifdef SPECT
    if (local_ele <= normPDF2((FLOAT_TWO+FLOAT_ONE+FLOAT_HALF)*CORstd, FLOAT_ZERO, CORstd)) {
#else
	if (local_ele <= THR) {
#endif
		return true;
	}
#endif
#if defined(BP) && defined(MASKBP) ///////////////////// APPLY BP MASK /////////////////////////
    if ((ii == 0) && (readMaskBP(maskBP, ind) == 0)) {
        return false;
    }
#endif ///////////////////// END BP MASK /////////////////////////
#if defined(SPECT) && defined(ATN)
	local_ele *= jelppi;
#endif
#if defined(FP) //////////////// FORWARD PROJECTION ////////////////
	denominator(ax, ind, local_ele, d_OSEM
#ifdef TOF //////////////// TOF ////////////////
				, element, TOFSum, DD, sigma_x, D, TOFWeights
#ifdef LISTMODE
        , TOFid
#endif
#endif //////////////// END TOF ////////////////
#ifdef N_RAYS //////////////// MULTIRAY
        , lor
#endif //////////////// END MULTIRAY
	);
#endif //////////////// END FORWARD PROJECTION ////////////////
#if defined(BP) //////////////// BACKWARD PROJECTION ////////////////
	rhs(local_ele * temp, ax, local_ind, d_output, no_norm, d_Summ
#ifdef TOF
				, element, sigma_x, D, DD, TOFSum, TOFWeights
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
	, const FLOAT element, const FLOAT sigma_x, float* D, const FLOAT DD, const FLOAT TOFSum, float* TOFWeights
#ifdef LISTMODE
	, const int TOFid
#endif
#endif
#if defined(SPECT) && defined(ATN)
	, FLOAT jelppi
#endif
#if defined(MASKBP) && defined(BP)
	, const int ii, MASKBPTYPE maskBP
#endif
#ifdef SPECT
    , const FLOAT coneOfResponseStdCoeffA, const FLOAT coneOfResponseStdCoeffB, const FLOAT coneOfResponseStdCoeffC, const FLOAT orth_ray_length_inv_signed
#endif
#ifdef N_RAYS
    , int lor
#endif
) {
	int uu = 0;
	bool breikki = false;
    FLOAT3 l; // Precomputed cross product elements
	const FLOAT v0 = center.x - s.x;
    l.z = diff.y * v0;
	const FLOAT apu1 = diff.z * v0;
	const int maksimiZ = CINT(d_Nz);
	const int minimiZ = 0;
	const int maksimiXY = CINT(d_N1);
	const int minimiXY = 0;
	int uu1 = 0, uu2 = 0;
#ifdef CRYSTZ
	int zz = tempk;
	for (zz = tempk; zz < maksimiZ; zz++) {
#else
	int zz = tempk;
#endif
		center.z = bz + CFLOAT(zz) * dz + dz / FLOAT_TWO;
		const FLOAT z0 = center.z - s.z;
		l.x = diff.x * z0 - apu1;
		l.y = diff.y * z0;
#ifdef CRYSTXY
		for (uu1 = temp2; uu1 < maksimiXY; uu1++) {
#else
		uu1 = temp2;
#endif
			center.y = b2 + CFLOAT(uu1) * d2 + d2 / FLOAT_TWO;
			breikki = orthogonalHelper3D(tempi, uu1, d_N2, d_N3, d_Nxy, zz, s, l, diff, orth_ray_length, center, bmin, bmax, Vmax, V, XY, ax, temp, 
#if defined(FP)
				d_OSEM
#else
				Summ, d_rhs_OSEM, no_norm
#endif
#ifdef TOF
				, element, sigma_x, D, DD, TOFSum, TOFWeights
#ifdef LISTMODE
				, TOFid
#endif
#endif
#if defined(SPECT) && defined(ATN)
				, jelppi
#endif
#if defined(MASKBP) && defined(BP)
				, ii, maskBP
#endif
#ifdef SPECT
                , coneOfResponseStdCoeffA, coneOfResponseStdCoeffB, coneOfResponseStdCoeffC, orth_ray_length_inv_signed
#endif
#ifdef N_RAYS
                    , lor
#endif
			);
#ifdef CRYSTXY
			if (breikki) {
				break;
			}
			uu++;
		}
#endif
#ifdef CRYSTXY
		for (uu2 = temp2 - 1; uu2 >= minimiXY; uu2--) {
			center.y = b2 + CFLOAT(uu2) * d2 + d2 / FLOAT_TWO;
			breikki = orthogonalHelper3D(tempi, uu2, d_N2, d_N3, d_Nxy, zz, s, l, diff, orth_ray_length, center, bmin, bmax, Vmax, V, XY, ax, temp, 
#if defined(FP)
				d_OSEM
#else
				Summ, d_rhs_OSEM, no_norm
#endif
#ifdef TOF
				, element, sigma_x, D, DD, TOFSum, TOFWeights
#ifdef LISTMODE
				, TOFid
#endif
#endif
#if defined(SPECT) && defined(ATN)
				, jelppi
#endif
#if defined(MASKBP) && defined(BP)
				, ii, maskBP
#endif
#ifdef SPECT
                , coneOfResponseStdCoeffA, coneOfResponseStdCoeffB, coneOfResponseStdCoeffC, orth_ray_length_inv_signed
#endif
#ifdef N_RAYS
                    , lor
#endif
			);
			if (breikki) {
				break;
			}
			uu++;
		}
#else
	uu2 = temp2 - 1;
#endif
#ifdef CRYSTZ
	if (uu1 == temp2 && uu2 == temp2 - 1 && breikki)
		break;
	}
	*k = zz - 1;
	for (zz = tempk - 1; zz >= minimiZ; zz--) {
		center.z = bz + CFLOAT(zz) * dz + dz / FLOAT_TWO;
		const FLOAT z0 = center.z - s.z;
		l.x = diff.x * z0 - apu1;
		l.y = diff.y * z0;
#ifdef CRYSTXY
		for (uu1 = temp2; uu1 < maksimiXY; uu1++) {
#else
		uu1 = temp2;
#endif
			center.y = b2 + CFLOAT(uu1) * d2 + d2 / FLOAT_TWO;
			breikki = orthogonalHelper3D(tempi, uu1, d_N2, d_N3, d_Nxy, zz, s, l, diff, orth_ray_length, center, bmin, bmax, Vmax, V, XY, ax, temp, 
#if defined(FP)
				d_OSEM
#else
				Summ, d_rhs_OSEM, no_norm
#endif
#ifdef TOF
				, element, sigma_x, D, DD, TOFSum, TOFWeights
#ifdef LISTMODE
				, TOFid
#endif
#endif
#if defined(SPECT) && defined(ATN)
				, jelppi
#endif
#if defined(MASKBP) && defined(BP)
				, ii, maskBP
#endif
#ifdef SPECT
                , coneOfResponseStdCoeffA, coneOfResponseStdCoeffB, coneOfResponseStdCoeffC, orth_ray_length_inv_signed
#endif
#ifdef N_RAYS
                    , lor
#endif
			);
#ifdef CRYSTXY
			if (breikki) {
				break;
			}
			uu++;
		}
#endif
#ifdef CRYSTXY
		for (uu2 = temp2 - 1; uu2 >= minimiXY; uu2--) {
			center.y = b2 + CFLOAT(uu2) * d2 + d2 / FLOAT_TWO;
			breikki = orthogonalHelper3D(tempi, uu2, d_N2, d_N3, d_Nxy, zz, s, l, diff, orth_ray_length, center, bmin, bmax, Vmax, V, XY, ax, temp, 
#if defined(FP)
				d_OSEM
#else
				Summ, d_rhs_OSEM, no_norm
#endif
#ifdef TOF
				, element, sigma_x, D, DD, TOFSum, TOFWeights
#ifdef LISTMODE
				, TOFid
#endif
#endif
#if defined(SPECT) && defined(ATN)
				, jelppi
#endif
#if defined(MASKBP) && defined(BP)
				, ii, maskBP
#endif
#ifdef SPECT
                , coneOfResponseStdCoeffA, coneOfResponseStdCoeffB, coneOfResponseStdCoeffC, orth_ray_length_inv_signed
#endif
#ifdef N_RAYS
                    , lor
#endif
			);
			if (breikki) {
				break;
			}
			uu++;
		}
#else
	uu2 = temp2 - 1;
#endif
		if (uu1 == temp2 && uu2 == temp2 - 1 && breikki)
			break;
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