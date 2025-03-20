
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
DEVICE float compute_element_orth_3D(const float xs, const float ys, const float zs, const float xl, const float yl, const float zl, const float crystal_size_z,
	const float xp) {
	const float x0 = xp - xs;

	// Cross product
#ifdef USEMAD
	const float y1 = FMAD(zl, x0, - xl);
	const float z1 = FMAD(-yl, x0, ys);
#else
	const float y1 = zl * x0 - xl;
	const float z1 = -yl * x0 + ys;
#endif
	const float normi = length(CMFLOAT3(zs, y1, z1));
#ifdef VOL
	return (normi / crystal_size_z);
#else
	return (1.f - normi / crystal_size_z);
#endif
}

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
) {
#if (defined(FP) || (defined(MASKBP) && defined(BP))) && defined(USEIMAGES)
	int3 ind;
		if (XY)
			ind = CMINT3(tempi, uu, zz);
		else
			ind = CMINT3(uu, tempi, zz);
#endif
	float local_ele = compute_element_orth_3D(s2, l3, l1, l2, diff1, diffZ, kerroin, center2);
#ifdef VOL
	if (local_ele >= bmax) {
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
 
// Both the orthogonal and volumme of intersection ray tracers loop through all the neighboring voxels of the current voxel
// Both also proceed through each X or Y slice, depending on the incident direction
// This function simply loops through each X or Y voxel and Z voxels in the current slice
// Forward or backward projection is computed in the helper function
DEVICE int orthDistance3D(const int tempi, const float diff1, const float diff2, const float diffZ, 
	const float center1, const float b2, const float d2, const float bz, const float dz, const float temp, int temp2, const int tempk, 
	const float s1, const float s2, const float sZ, const uint d_Nxy, const float kerroin,
	const uint d_N1, const uint d_N2, const uint d_N3, const uint d_Nz, const float bmin, 
	const float bmax, const float Vmax, CONSTANT float* V, const bool XY, float* ax, 
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
	for (int zz = tempk; zz < maksimiZ; zz++) {
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
	for (int zz = tempk - 1; zz >= minimiZ; zz--) {
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
#endif
	return uu;
}