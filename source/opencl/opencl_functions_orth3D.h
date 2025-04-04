
/*******************************************************************************************************************************************
* Special functions for the 3D orthogonal distance-based ray tracer and volume of intersection ray tracer.
*
* Copyright (C) 2020-2024 Ville-Veikko Wettenhovi
*
* This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************************************************************************/

// Compute the orthogonal distance from the ray to the current voxel (center)
// For orthogonal distance-based ray tracer, the distance is normalized
DEVICE float compute_element_orth_3D(const float xs, const float ys, const float zs, const float xl, const float yl, const float zl, const float crystal_size_z, const float xp) {
	const float x0 = xp - xs;

	// Cross product
#ifdef USEMAD
	const float y1 = FMAD(zl, x0, - xl);
	const float z1 = FMAD(-yl, x0, ys);
#else
	const float y1 = zl * x0 - xl;
	const float z1 = -yl * x0 + ys;
#endif
	const float norm1 = length(CMFLOAT3(zs, y1, z1));
#ifdef SPECT
    const float norm2 = length(CMFLOAT3(x0, yl, zl));
    const float d = norm1 / norm2;
    return d;
#else
#ifdef VOL
	return (norm1 / crystal_size_z);
#else
	return (1.f - (norm1 / crystal_size_z));
#endif
#endif
}

#ifdef SPECT
DEVICE float compute_element_parallel_3D(const float v0x, const float v0y, const float v0z, const float v1x, const float v1y, const float v1z, const float px, const float py, const float pz) {
    // In this function the ray is defined as v0+t*v1, where v0 is the source end of the ray and v1x for example is detectors.xd-detectors.xs
#ifdef USEMAD
    const float dot1 = FMAD(v1x, (px-v0x), FMAD(v1y, (py-v0y), FMAD(v1z, (pz-v0z), 0.f)));
    const float dot2 = FMAD(v1x, v1x, FMAD(v1y, v1y, FMAD(v1z, v1z, 0.f)));
    const float t = FMAD(-1.f, DIVIDE(dot1, dot2), 1.f);
    const float rayLength = FMAD(t, length(CMFLOAT3(v1x, v1y, v1z)), 0.f);
#else
    const float dot1 = v1x*(px-v0x)+v1y*(py-v0y)+v1z*(pz-v0z); // v1 * (p-v0)
    const float dot2 = v1x*v1x+v1y*v1y+v1z*v1z; // v1 * v1
    const float t = 1.f - dot1 / dot2; // 1-t as SPECT collimator response is measured from the collimator-detector interface
    const float rayLength = t * length(CMFLOAT3(v1x, v1y, v1z));
#endif
    return (rayLength);
}

#define _2PI 0.3989423f
DEVICE float normPDF2(const float x, const float mu, const float sigma) {
#ifdef USEMAD
    const float a = FMAD(-1.f, DIVIDE(mu, sigma), DIVIDE(x, sigma));
    return POWR(DIVIDE(_2PI, sigma), 2.f) * EXP(-0.5f * a * a);
#else
	const float a = (x - mu) / sigma;
	return _2PI * _2PI / sigma / sigma * EXP(-0.5f * a * a);
#endif
}
#endif

// compute voxel index, orthogonal distance based or volume of intersection ray tracer
DEVICE LONG compute_ind_orth_3D(const uint tempi, const uint tempijk, const int tempk, const uint d_N, const uint Nyx) {
	LONG local_ind = CLONG_rtz(tempi * d_N + tempijk) + CLONG_rtz(tempk) * CLONG_rtz(Nyx);
	return local_ind;
}

// This function computes either the forward or backward projection for the current voxel
// The normalized orthogonal distance or the volume of the (spherical) voxel is computed before the forward or backward projections
DEVICE bool orthogonalHelper3D(const int tempi, const int uu, const uint d_N2, const uint d_N3, const uint d_Nxy, const int zz, const float s2, const float l3, const float l1, const float l2,
	const float diff1, const float diffZ, const float kerroin, const float center2, const float bmin, const float bmax, const float Vmax,
	CONSTANT float* V, const bool XY, float* ax, const float temp, 
#if defined(FP)
	IMTYPE d_OSEM
#else
	CLGLOBAL CAST* d_Summ, CLGLOBAL CAST* d_output, const bool no_norm
#endif
#ifdef TOF
	, const float element, const float sigma_x, float* D, const float DD, CONSTANT float* TOFCenter, const float TOFSum
#ifdef LISTMODE
	, const int TOFid
#endif
#endif
#if defined(SPECT) && defined(ATN)
	, float jelppi
#endif
#if defined(MASKBP) && defined(BP)
#ifdef USEIMAGES
#ifdef MASKBP3D
	, const int ii, IMAGE3D maskBP
#else
	, const int ii, IMAGE2D maskBP
#endif
#else
, const int ii, const CLGLOBAL uchar* CLRESTRICT maskBP
#endif
#endif
#ifdef SPECT
    , const float coneOfResponseStdCoeffA, const float coneOfResponseStdCoeffB, const float coneOfResponseStdCoeffC, const float2 crXY, const float s1, const float sZ, const float diff2, const float center1, const float centerZ
#endif
) {
#if (defined(FP) || (defined(MASKBP) && defined(BP))) && defined(USEIMAGES)
	int3 ind;
		if (XY)
			ind = CMINT3(tempi, uu, zz);
		else
			ind = CMINT3(uu, tempi, zz);
#endif
#ifdef SPECT
    float d_orth = compute_element_orth_3D(s2, l3, l1, l2, diff1, diffZ, kerroin, center2);
    float d_parallel = compute_element_parallel_3D(s1, s2, sZ, diff1, diff2, diffZ, center1, center2, centerZ);
    if (d_parallel < 0) { // Voxel behind detector
        return true;
    }
#ifdef USEMAD
    float CORstd = DIVIDE(SQRT(FMAD(1.f, POWR(FMAD(coneOfResponseStdCoeffA, d_parallel, coneOfResponseStdCoeffB), 2.f), POWR(coneOfResponseStdCoeffC, 2.f))), (2.f*SQRT(2.f*LOG(2.f))));
#else
    float CORstd = sqrt(pow(coneOfResponseStdCoeffA*d_parallel+coneOfResponseStdCoeffB, 2.f)+pow(coneOfResponseStdCoeffC, 2.f)) / (2.f*sqrt(2.f*log(2.f)));
#endif
    float local_ele = normPDF2(d_orth, 0.f, CORstd);
#else
	float local_ele = compute_element_orth_3D(s2, l3, l1, l2, diff1, diffZ, kerroin, center2);
#endif
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
    if (local_ele <= normPDF2(3.5*CORstd, 0.f, CORstd)) {
#else
	if (local_ele <= THR) {
#endif
		return true;
	}
#endif
#if defined(SPECT) && defined(ATN)
	local_ele *= jelppi;
#endif
#ifdef BP
#if defined(MASKBP)
	if (ii == 0) {
#ifdef USEIMAGES
#ifdef CUDA
#ifdef MASKBP3D
		const int maskVal = tex3D<unsigned char>(maskBP, ind.x, ind.y, ind.z);
#else
		const int maskVal = tex2D<unsigned char>(maskBP, ind.x, ind.y);
#endif
#else
#ifdef MASKBP3D
		const int maskVal = read_imageui(maskBP, sampler_MASK, (int4)(ind.x, ind.y, ind.z, 0)).w;
#else
		const int maskVal = read_imageui(maskBP, sampler_MASK, (int2)(ind.x, ind.y)).w;
#endif
#endif
#else
#ifdef MASKBP3D
		const int maskVal = maskBP[tempi * d_N2 + uu * d_N3 + zz * d_Nxy];
#else
		const int maskVal = maskBP[tempi * d_N2 + uu * d_N3];
#endif
#endif
		if (maskVal == 0)
			return false;
	}
#endif
#endif
#if defined(BP) || (defined(FP) && !defined(USEIMAGES))
	LONG local_ind = 0;
	local_ind = compute_ind_orth_3D(CUINT(tempi), uu * d_N3, (zz), d_N2, d_Nxy);
#if defined(FP) && !defined(USEIMAGES)
	const LONG ind = CLONG_rtz(local_ind);
#endif
#endif
#if defined(FP) //////////////// FORWARD PROJECTION ////////////////
	denominator(ax, ind, local_ele, d_OSEM
#ifdef TOF //////////////// TOF ////////////////
				, element, TOFSum, DD, TOFCenter, sigma_x, D
#ifdef LISTMODE
				, TOFid
#endif
#endif //////////////// END TOF ////////////////
	);
#endif
#if defined(BP) //////////////// BACKWARD PROJECTION ////////////////
	rhs(local_ele * temp, ax, local_ind, d_output, no_norm, d_Summ
#ifdef TOF
				, element, sigma_x, D, DD, TOFCenter, TOFSum
#ifdef LISTMODE
				, TOFid
#endif
#endif
				);
#endif
	return false;
}
 
// Both the orthogonal and volume of intersection ray tracers loop through all the neighboring voxels of the current voxel
// Both also proceed through each X or Y slice, depending on the incident direction
// This function simply loops through each X or Y voxel and Z voxels in the current slice
// Forward or backward projection is computed in the helper function
DEVICE int orthDistance3D(const int tempi, const float diff1, const float diff2, const float diffZ, 
	const float center1, const float b2, const float d2, const float bz, const float dz, const float temp, int temp2, const int tempk, 
	const float s1, const float s2, const float sZ, const uint d_Nxy, const float kerroin,
	const uint d_N1, const uint d_N2, const uint d_N3, const uint d_Nz, const float bmin, 
	const float bmax, const float Vmax, CONSTANT float* V, const bool XY, float* ax, const bool preStep, int* k, const int ku, 
#if defined(FP)
	IMTYPE d_OSEM
#else
	const bool no_norm, CLGLOBAL CAST* Summ, CLGLOBAL CAST* d_rhs_OSEM 
#endif
#ifdef TOF
	, const float element, const float sigma_x, float* D, const float DD, CONSTANT float* TOFCenter, const float TOFSum
#ifdef LISTMODE
	, const int TOFid
#endif
#endif
#if defined(SPECT) && defined(ATN)
	, float jelppi
#endif
#if defined(MASKBP) && defined(BP)
#ifdef USEIMAGES
#ifdef MASKBP3D
	, const int ii, IMAGE3D maskBP
#else
	, const int ii, IMAGE2D maskBP
#endif
#else
	, const int ii, const CLGLOBAL uchar* CLRESTRICT maskBP
#endif
#endif
#ifdef SPECT
    , const float coneOfResponseStdCoeffA, const float coneOfResponseStdCoeffB, const float coneOfResponseStdCoeffC, const float2 crXY
#endif
) {
	int uu = 0;
	bool breikki = false;
	// y0
	const float v0 = center1 - s1;
	// xl * y0
	const float l3 = diff1 * v0;
	// zl * y0
	const float apu1 = diffZ * v0;
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
		const float centerZ = bz + CFLOAT(zz) * dz + dz / 2.f;
		const float z0 = centerZ - sZ;
		// x1 = yl * z0 - zl * y0
		const float l1 = diff2 * z0 - apu1;
		// xl * z0
		const float l2 = diff1 * z0;
#ifdef CRYSTXY
		for (uu1 = temp2; uu1 < maksimiXY; uu1++) {
#else
		uu1 = temp2;
#endif
			const float center2 = b2 + CFLOAT(uu1) * d2 + d2 / 2.f;
			breikki = orthogonalHelper3D(tempi, uu1, d_N2, d_N3, d_Nxy, zz, s2, l3, l1, l2, diff2, diffZ, kerroin, center2, bmin, bmax, Vmax, V, 
				XY, ax, temp, 
#if defined(FP)
				d_OSEM
#else
				Summ, d_rhs_OSEM, no_norm
#endif
#ifdef TOF
				, element, sigma_x, D, DD, TOFCenter, TOFSum
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
                , coneOfResponseStdCoeffA, coneOfResponseStdCoeffB, coneOfResponseStdCoeffC, crXY, s1, sZ, diff1, center1, centerZ
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
			const float center2 = b2 + CFLOAT(uu2) * d2 + d2 / 2.f;
			breikki = orthogonalHelper3D(tempi, uu2, d_N2, d_N3, d_Nxy, zz, s2, l3, l1, l2, diff2, diffZ, kerroin, center2, bmin, bmax, Vmax, V, 
				XY, ax, temp, 
#if defined(FP)
				d_OSEM
#else
				Summ, d_rhs_OSEM, no_norm
#endif
#ifdef TOF
				, element, sigma_x, D, DD, TOFCenter, TOFSum
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
                , coneOfResponseStdCoeffA, coneOfResponseStdCoeffB, coneOfResponseStdCoeffC, crXY, s1, sZ, diff1, center1, centerZ
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
		const float centerZ = bz + CFLOAT(zz) * dz + dz / 2.f;
		const float z0 = centerZ - sZ;
		const float l1 = diff2 * z0 - apu1;
		const float l2 = diff1 * z0;
#ifdef CRYSTXY
		for (uu1 = temp2; uu1 < maksimiXY; uu1++) {
#else
		uu1 = temp2;
#endif
			const float center2 = b2 + CFLOAT(uu1) * d2 + d2 / 2.f;
			breikki = orthogonalHelper3D(tempi, uu1, d_N2, d_N3, d_Nxy, zz, s2, l3, l1, l2, diff2, diffZ, kerroin, center2, bmin, bmax, Vmax, V, 
				XY, ax, temp, 
#if defined(FP)
				d_OSEM
#else
				Summ, d_rhs_OSEM, no_norm
#endif
#ifdef TOF
				, element, sigma_x, D, DD, TOFCenter, TOFSum
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
                , coneOfResponseStdCoeffA, coneOfResponseStdCoeffB, coneOfResponseStdCoeffC, crXY, s1, sZ, diff1, center1, centerZ
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
			const float center2 = b2 + CFLOAT(uu2) * d2 + d2 / 2.f;
			breikki = orthogonalHelper3D(tempi, uu2, d_N2, d_N3, d_Nxy, zz, s2, l3, l1, l2, diff2, diffZ, kerroin, center2, bmin, bmax, Vmax, V, 
				XY, ax, temp, 
#if defined(FP)
				d_OSEM
#else
				Summ, d_rhs_OSEM, no_norm
#endif
#ifdef TOF
				, element, sigma_x, D, DD, TOFCenter, TOFSum
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
                , coneOfResponseStdCoeffA, coneOfResponseStdCoeffB, coneOfResponseStdCoeffC, crXY, s1, sZ, diff1, center1, centerZ
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