
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

#ifdef SPECT // Threshold for ODRT collimator response / tube width cutoff
#define THR normPDF2(3.5*CORstd, 0.f, CORstd)
#else
#define THR 0.01f
#endif

#ifndef PTYPE4
// Compute the orthogonal distance from the ray to the current voxel (center)
// For orthogonal distance-based ray tracer, the distance is normalized
// See for example https://math.stackexchange.com/questions/2353288/point-to-line-distance-in-3d-using-cross-product/2355960#2355960
DEVICE float compute_element_orth_3D(
    const float3 s, // Source position
    const float3 l, // Precomputed cross product elements
    const float3 diff, // Spans the ray
    const float3 center, // Center of the voxel
    const float orth_ray_length // Precomputed length |diff|, in PET scaled with tube width
) {
    // In this function the ray is defined as v0+t*v1, where v0 is the source end of the ray and v0+v1 is the detector end of the ray. Thus 0<=t<=1. The orthogonal distance d from the voxel center p to the ray is equal to d = |v1 x (v0-p)| / |v1|. The term |v1| is precomputed, and in PET reconstruction includes the crystal size.

    // Precomputed:
    //l.x = diff.x * (center.z - s.z) - diff.z  * (center.x - s.x);
    //l.y = diff.y * (center.z - s.z);
    //l.z = diff.y * (center.x - s.x);

#ifdef USEMAD
    const float x0 = FMAD(-1.f, s.y, center.y);
	const float y1 = FMAD(diff.z, x0, - l.y);
	const float z1 = FMAD(-diff.x, x0, l.z);
#else
    const float x0 = center.y - s.y;
	const float y1 = diff.z * x0 - l.y;
	const float z1 = -diff.x * x0 + l.z;
#endif
	const float norm1 = length(CMFLOAT3(l.x, y1, z1));
    return 
#if !defined(VOL) && !defined(SPECT)
    1.f - 
#endif
    norm1 / orth_ray_length;
}

#ifdef SPECT
// In this function the ray is defined as v0+t*v1, where v0 is the source end of the ray and v0+v1 is the detector end of the ray. Thus 0<=t<=1. The parallel distance d is computed as the distance from the source to the projection of point p onto the ray. That is, d = v1*(p-v0) / (v1*v1) * |v1|
//return dot(v1, p - v0) * length(v1) / dot(v1, v1);
DEVICE float compute_element_parallel_3D(const float3 v0, const float3 v1, const float3 p, const float orth_ray_length_inv_signed) {
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
    const float3 s,
    const float3 l,
	const float3 diff, 
    const float orth_ray_length, const float3 center, const float bmin, const float bmax, const float Vmax, CONSTANT float* V, const bool XY, float* ax, const float temp, 
#if defined(FP)
	IMTYPE d_OSEM
#else
	CLGLOBAL CAST* d_Summ, CLGLOBAL CAST* d_output, const bool no_norm
#endif
#ifdef TOF
	, const float element, const float sigma_x, float* D, const float DD, const float TOFSum, float* TOFWeights
#ifdef LISTMODE
	, const int TOFid
#endif
#endif
#if defined(SPECT) && defined(ATN)
	, float jelppi
#endif
#if defined(MASKBP) && defined(BP)
    , const int ii, MASKBPTYPE maskBP
#endif
#ifdef SPECT
    , const float coneOfResponseStdCoeffA, const float coneOfResponseStdCoeffB, const float coneOfResponseStdCoeffC, const float orth_ray_length_inv_signed
#endif
#ifdef N_RAYS
    , int lor
#endif
) {
#ifdef SPECT // Check for voxel behind detector
    float d_parallel = compute_element_parallel_3D(s, diff, center, orth_ray_length_inv_signed);
    if (d_parallel < 0) { 
        return false;
    }
#endif
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
    float local_ele = compute_element_orth_3D(s, l, diff, center, orth_ray_length);
#ifdef SPECT ////////////////////////// SPECT ////////////////////
    float t = FMAD(coneOfResponseStdCoeffA, d_parallel, coneOfResponseStdCoeffB);   // A*d_parallel + B
    float t2 = t * t;
    float c2 = coneOfResponseStdCoeffC * coneOfResponseStdCoeffC;
    float CORstd = SQRT(t2 + c2) * INV_2SQRT2LN2;
    local_ele = normPDF2(local_ele, 0.f, CORstd);
#endif //////////// END SPECT / NOT SPECT /////////////
#ifdef VOL
	if (local_ele > bmax) {
		return true;
	}
	if (local_ele < bmin)
		local_ele = Vmax;
	else
		local_ele = V[CUINT_rte((local_ele - bmin) * CC)];
#else
	if (local_ele <= THR) {
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
    const float3 diff, // detectors.[x/y/z]d - detectors.[x/y/z]s
    float3 center, // center of the current voxel
    const float3 s, // detectors.[x/y/z]s
    const float b2, const float d2, const float bz, const float dz, const float temp, int temp2, const int tempk, 
	const uint d_Nxy, const float orth_ray_length,
	const uint d_N1, const uint d_N2, const uint d_N3, const uint d_Nz, const float bmin, 
	const float bmax, const float Vmax, CONSTANT float* V, const bool XY, float* ax, const bool preStep, int* k, const int ku, 
#if defined(FP)
	IMTYPE d_OSEM
#else
	const bool no_norm, CLGLOBAL CAST* Summ, CLGLOBAL CAST* d_rhs_OSEM 
#endif
#ifdef TOF
	, const float element, const float sigma_x, float* D, const float DD, const float TOFSum, float* TOFWeights
#ifdef LISTMODE
	, const int TOFid
#endif
#endif
#if defined(SPECT) && defined(ATN)
	, float jelppi
#endif
#if defined(MASKBP) && defined(BP)
	, const int ii, MASKBPTYPE maskBP
#endif
#ifdef SPECT
    , const float coneOfResponseStdCoeffA, const float coneOfResponseStdCoeffB, const float coneOfResponseStdCoeffC, const float orth_ray_length_inv_signed
#endif
#ifdef N_RAYS
    , int lor
#endif
) {
	int uu = 0;
	bool breikki = false;
    float3 l; // Precomputed cross product elements
	// y0
	const float v0 = center.x - s.x;
	// xl * y0
    l.z = diff.y * v0;
	// zl * y0
	const float apu1 = diff.z * v0;
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
		// z0
		center.z = bz + CFLOAT(zz) * dz + dz / 2.f;
		const float z0 = center.z - s.z;
		// x1 = yl * z0 - zl * y0
		l.x = diff.x * z0 - apu1;
		// xl * z0
		l.y = diff.y * z0;
#ifdef CRYSTXY
		for (uu1 = temp2; uu1 < maksimiXY; uu1++) {
#else
		uu1 = temp2;
#endif
			center.y = b2 + CFLOAT(uu1) * d2 + d2 / 2.f;
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
			center.y = b2 + CFLOAT(uu2) * d2 + d2 / 2.f;
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
		center.z = bz + CFLOAT(zz) * dz + dz / 2.f;
		const float z0 = center.z - s.z;
		l.x = diff.x * z0 - apu1;
		l.y = diff.y * z0;
#ifdef CRYSTXY
		for (uu1 = temp2; uu1 < maksimiXY; uu1++) {
#else
		uu1 = temp2;
#endif
			center.y = b2 + CFLOAT(uu1) * d2 + d2 / 2.f;
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
			center.y = b2 + CFLOAT(uu2) * d2 + d2 / 2.f;
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
			*k = max(*k, zz) - 1;
		else
			*k = min(*k, zz) + 1;
	}
#endif
	return uu;
}

#endif