/**************************************************************************
* Header file for the various functions required by the improved Siddon
* ray tracer.
*
* This file can be used to compute the forward and/or backward projections
* in any C++-code. Use the below function projectorType123Implementation4
* to compute either the forward or backward projection. Use paramStruct
* struct to input the necessary and optinal input parameters.
*
* Copyright (C) 2019-2024 Ville-Veikko Wettenhovi, Niilo Saarlemo
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <https://www.gnu.org/licenses/>.
***************************************************************************/
#pragma once

#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <numeric>
#include <time.h>
#include <thread>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#ifndef nChunks
#define nChunks 1000
#endif
#endif
#ifdef MATLAB
#include "mex.h"
#endif

// TOF constants, the same as in the GPU kernels (general_opencl_functions.h)
// Distance in sigmas after which the Gaussian is considered to be zero
#define TOF_CUTOFF 6.
#define TOF_INV_SQRT_2PI 0.3989422804014327

inline void setThreads() {
#ifdef _OPENMP
	if (omp_get_max_threads() == 1) {
		int n_threads = std::thread::hardware_concurrency();
		omp_set_num_threads(n_threads);
	}
	else
		omp_set_num_threads(omp_get_max_threads());
#endif
}

template <typename T> T sign(T val) {
	return (T(0) < val) - (val < T(0));
}

template <typename T>
struct Det {
	T xd = 0., xs = 0., yd = 0., ys = 0., zd = 0., zs = 0.;
};

template <typename T>
struct float2 {
	T x = 0.f, y = 0.f;
};

template <typename T>
struct float3 {
	T x = 0.f, y = 0.f, z = 0.f;
};

template <typename T>
void MFLOAT2(T x, T y, float2<T>& xy) {
	xy.x = x;
	xy.y = y;
}

template <typename T>
void MFLOAT3(T x, T y, T z, float3<T>& xyz) {
	xyz.x = x;
	xyz.y = y;
	xyz.z = z;
}

template <typename T>
float2<T> operator+(float2<T> a, float2<T> b) {
	float2<T> out;
	out.x = a.x + b.x;
	out.y = a.y + b.y;
	return out;
}

template <typename T>
float2<T> operator+(const T a, float2<T> b) {
	float2<T> out;
	out.x = a + b.x;
	out.y = a + b.y;
	return out;
}

template <typename T>
float2<T> operator+(float2<T> a, const T b) {
	float2<T> out;
	out.x = a.x + b;
	out.y = a.y + b;
	return out;
}

template <typename T>
float2<T> operator-(const T a, float2<T> b) {
	float2<T> out;
	out.x = a - b.x;
	out.y = a - b.y;
	return out;
}

template <typename T>
float2<T> operator*(const T a, float2 <T>b) {
	float2<T> out;
	out.x = a * b.x;
	out.y = a * b.y;
	return out;
}

template <typename T>
float2<T> operator-(float2<T> a, float2<T> b) {
	float2<T> out;
	out.x = a.x - b.x;
	out.y = a.y - b.y;
	return out;
}

template <typename T>
float2<T> operator*(float2<T> a, float2<T> b) {
	float2<T> out;
	out.x = a.x * b.x;
	out.y = a.y * b.y;
	return out;
}

template <typename T>
float2<T> operator/(float2<T> a, float2<T> b) {
	float2<T> out;
	out.x = a.x / b.x;
	out.y = a.y / b.y;
	return out;
}

template <typename T>
float2<T> operator/(float2<T> a, const float b) {
	float2<T> out;
	out.x = a.x / b;
	out.y = a.y / b;
	return out;
}

template <typename T>
float2<T> fabs2(float2<T> a) {
	float2<T> out;
	out.x = std::fabs(a.x);
	out.y = std::fabs(a.y);
	return out;
}

template <typename T>
float3<T> operator-(float3<T> a, const float b) {
	float3<T> out;
	out.x = a.x - b;
	out.y = a.y - b;
	out.z = a.z - b;
	return out;
}

template <typename T>
float3<T> operator-(float3<T> a) {
	float3<T> out;
	out.x = -a.x;
	out.y = -a.y;
	out.z = -a.z;
	return out;
}

template <typename T>
float3<T> operator*(float3<T> a, float3<T> b) {
	float3<T> out;
	out.x = a.x * b.x;
	out.y = a.y * b.y;
	out.z = a.z * b.z;
	return out;
}

template <typename T>
float3<T> operator/(float3<T> a, const float b) {
	float3<T> out;
	out.x = a.x / b;
	out.y = a.y / b;
	out.z = a.z / b;
	return out;
}

template <typename T>
float3<T> exp3(float3<T> a) {
	float3<T> out;
	out.x = std::exp(a.x);
	out.y = std::exp(a.y);
	out.z = std::exp(a.z);
	return out;
}

template <typename T>
struct paramStruct {
	// Subset type
	int subsetType = 1;
	// Number of subsets
	uint32_t subsets = 1;
	// Current subset number (zero-based indexing)
	uint32_t currentSubset = 0;
	// Horizontal number of detector elements (REQUIRED)
	uint32_t size_y = 1;
	// Vertical number of detector elements (REQUIRED)
	uint32_t size_x = 1;
	// Voxel width (x-direction) (REQUIRED)
	T dx = static_cast<T>(1.);
	// Voxel width (y-direction) (REQUIRED)
	T dy = static_cast<T>(1.);
	// Voxel width (z-direction) (REQUIRED)
	T dz = static_cast<T>(1.);
	// Number of voxels in the x-direction (REQUIRED)
	uint32_t Nx = 1U;
	// Number of voxels in the y-direction (REQUIRED)
	uint32_t Ny = 1U;
	// Number of voxels in the z-direction (REQUIRED)
	uint32_t Nz = 1U;
	// Coordinate of the (lower left) edge of the FOV with respect the origin (x-direction) (REQUIRED)
	T bx = static_cast<T>(1.);
	// Coordinate of the (lower left) edge of the FOV with respect the origin (y-direction) (REQUIRED)
	T by = static_cast<T>(1.);
	// Coordinate of the (lower left) edge of the FOV with respect the origin (z-direction) (REQUIRED)
	T bz = static_cast<T>(1.);
	// If true, the forward projection mask is used
	bool useMaskFP = false;
	// Forward projection mask
	uint8_t* maskFP = nullptr;
	// Detectors per one PET ring (required for listmode sensitivity image)
	uint32_t det_per_ring = 1;
	// number of PET rings (required for listmode sensitivity image)
	int32_t rings = 1;
	// Number of transaxial rays (multi-ray reconstruction only)
	uint16_t nRays2D = 1;
	// Number of axial rays (multi-ray reconstruction only)
	uint16_t nRays3D = 1;
	// Whether list-mode data is used (1 for list-mode, 0 for sinogram)
	uint8_t listMode = 0;
	// If true, normalization correction is used
	bool normalizationCorrection = false;
	// If true, attenuation correction is used
	bool attenuationCorrection = false;
	// If true, attenuation images (e.g. scaled CT images) are used
	// If false, measurement (LOR) based attenuation correction is used (e.g. blank/transmission scan)
	bool CTAttenuation = true;
	// If true, multiplicative scatter correction is used  (i.e. the scatter data is multiplied with the system matrix)
	bool scatterCorrectionMult = false;
	// Normalization coefficients
	T* normCoefs = nullptr;
	// Scatter coefficients
	T* scatterCoefs = nullptr;
	// Attenuation coefficients
	T* atten = nullptr;
	// Detector pitch (size) in the transaxial (XY) direction (REQUIRED FOR MULTI-RAY OR CT)
	T dPitchXY = static_cast<T>(1.);
	// Detector pitch (size) in the saxial (Z) direction (REQUIRED FOR MULTI-RAY OR CT)
	T dPitchZ = static_cast<T>(1.);
	// XY-detector indices for subset types 3, 6 and 7, i.e. the index for the transaxial detector coordinate for the current LOR
	uint32_t* xy_index = nullptr;
	// Same as above, but for axial direction
	uint16_t* z_index = nullptr;
	// Projector type, only 1 (improved Siddon) is currently supported
	int projType = 1;
	// (CT ONLY) if the detector panel is not exactly centered, set this to true
	bool pitch = false;
	// (PET ONLY) is raw data used
	bool raw = false;
	// If set as true, the sensitivity image will NOT be computed
	bool noSensImage = true;
	// If true, TOF data is assumed
	bool TOF = false;
	// Number of TOF bins (this should always be at least 1)
	uint32_t nBins = 1;
	// TOF FWHM (REQUIRED FOR TOF)
	T sigma_x = static_cast<T>(1.);
	// A vector of TOF time distances from the center bin, i.e. 0, +1, -1, +2, -2, etc. (REQUIRED FOR TOF)
	T* TOFCenters = nullptr;
	// Number of projections/sinograms (in the current subset) (REQUIRED FOR CT OR PET WITH SUBSET TYPE >= 8)
	uint32_t nProjections = 1U;
	// (PET ONLY) Global correction factor
	T globalFactor = static_cast<T>(1.);
	// Small parameter to prevent division by zero
	T epps = static_cast<T>(1e-8);
	// Are multiple layers used in PET
	bool nLayers = false;
	// Number of PET detector layers (list-mode sensitivity image only, GPU NLAYERS). The detector coordinates x contain the detectors of
	// layer l at the offset l * det_per_ring and z the rings of layer l at the offset l * rings
	uint32_t nLayerCount = 1U;
	//
	uint64_t nMeas = 0;
	// if true, then computes the sensitivity image for listmode data
	bool computeSensIm = false;
	// If true, uses a unsigned 8-bit mask in backward projection
	bool useMaskBP = false;
	// Backward projection mask
	uint8_t* maskBP = nullptr;
	// SPECT ray shifts, detector end (collimator model)
	T* rayShiftsDetector = nullptr;
	// SPECT ray shifts, source end (collimator model)
	T* rayShiftsSource = nullptr;
    // SPECT cone of response standard deviation parameters
    T coneOfResponseStdCoeffA = 0.01;
    T coneOfResponseStdCoeffB = 0.01;
    T coneOfResponseStdCoeffC = 0.01;
	// (PET AND SPECT) If true, the rays are normalized with their total length (GPU -DTOTLENGTH), otherwise with the length of the ray inside the FOV.
	// For SPECT the total length is the length of the ray inside the support (see extendRayToSupport)
	bool useTotLength = true;
	// (SPECT ONLY) Detector head index of each projection (GPU: d_detectorVector), indexed with the same projection index as the detector
	// coordinates x. Selects the ray shifts of the head. If nullptr, every projection uses the head 0
	const uint32_t* detectorVector = nullptr;
	// (SPECT ONLY) The support (box or elliptic cylinder) the rays are clipped to (GPU: extendRayToEllipse). ellipsePower > 1e20 (e.g. Inf
	// or FLT_MAX) is a box, 2 an elliptic cylinder (with the axial extent ellipseCenterZ +- ellipseRadiusZ). Any other power leaves the
	// rays unclipped, as on the GPU. A non-positive radius (the default) uses the image volume as the box support
	T ellipseCenterX = static_cast<T>(0.);
	T ellipseCenterY = static_cast<T>(0.);
	T ellipseCenterZ = static_cast<T>(0.);
	T ellipseRadiusX = static_cast<T>(-1.);
	T ellipseRadiusY = static_cast<T>(-1.);
	T ellipseRadiusZ = static_cast<T>(-1.);
	T ellipsePower = static_cast<T>(1e30);
};

// Compute the Euclidean norm of a vector
template <typename T>
inline T norm(const T x, const T y, const T z) {
	return sqrt(x * x + y * y + z * z);
}

// Get the detector coordinates for the current list-mode measurement
template <typename T>
inline void get_detector_coordinates_raw(const T* x, const T* z, Det<T>& detectors,
	const size_t ll, const uint8_t list_mode_format, const uint16_t n_rays, const uint16_t n_rays3D, 
	const int ix, const int iy, const int izS, const int izD,
	const int lorZ, const int lorXY, const T crz, const T crxy, const bool Sens = false) {
	if (list_mode_format == 1 && !Sens) {
		detectors.xs = x[ll * 6];
		detectors.ys = x[ll * 6 + 1];
		detectors.zs = x[ll * 6 + 2];
		detectors.xd = x[ll * 6 + 3];
		detectors.yd = x[ll * 6 + 4];
		detectors.zd = x[ll * 6 + 5];
	}
	else {
		detectors.xs = x[ix * 2];
		detectors.ys = x[ix * 2 + 1];
		detectors.xd = x[iy * 2];
		detectors.yd = x[iy * 2 + 1];
		detectors.zs = z[izS];
		detectors.zd = z[izD];
	}
	if (n_rays3D > 1)
		multirayCoordinateShiftZ(detectors, lorZ, crz, n_rays3D);
	if (n_rays > 1)
		multirayCoordinateShiftXY(detectors, lorXY, crxy, n_rays);
}

// Get the detector coordinates for the current sinogram bin
template <typename T>
inline void get_detector_coordinates(const T* x, const T* z, const uint32_t size_x, const uint32_t size_y, Det<T>& detectors, const uint32_t* xy_index,
	const uint16_t* z_index, const size_t oo, const int subsetType, const uint32_t subsets, const int64_t ix, const int64_t iy, const int64_t iz,
	const uint16_t n_rays, const uint16_t n_rays3D, const int lorZ, const int lorXY, const T crz, const T crxy, const bool nLayers = false) {
	// Sinogram data
	if (subsetType >= 8 || subsets == 1 || subsetType == 1 || subsetType == 2 || subsetType == 4 || subsetType == 5 || subsetType == 0) {
		const int64_t id = (ix + iy * size_x) * 4LL;
		if (nLayers) {
			// With multiple layers z has three values per sinogram, the third one is the layer combination
			// (the transaxial coordinates are stored as one block of size_x * size_y bins per layer combination)
			const int64_t idz = iz * 3LL;
			const int64_t layerOffset = static_cast<int64_t>(z[idz + 2]) * static_cast<int64_t>(size_x) * static_cast<int64_t>(size_y) * 4LL;
			detectors.xs = x[id + layerOffset];
			detectors.ys = x[id + 1 + layerOffset];
			detectors.xd = x[id + 2 + layerOffset];
			detectors.yd = x[id + 3 + layerOffset];
			detectors.zs = z[idz];
			detectors.zd = z[idz + 1];
		}
		else {
			const int64_t idz = iz * 2LL;
			detectors.xs = x[id];
			detectors.ys = x[id + 1];
			detectors.xd = x[id + 2];
			detectors.yd = x[id + 3];
			detectors.zs = z[idz];
			detectors.zd = z[idz + 1];
		}
	}
	else {
		const uint32_t ind = xy_index[oo] * 4U;
		if (nLayers) {
			const uint32_t indz = (uint32_t)z_index[oo] * 3U;
			const int64_t layerOffset = static_cast<int64_t>(z[indz + 2]) * static_cast<int64_t>(size_x) * static_cast<int64_t>(size_y) * 4LL;
			detectors.xs = x[ind + layerOffset];
			detectors.ys = x[ind + 1 + layerOffset];
			detectors.xd = x[ind + 2 + layerOffset];
			detectors.yd = x[ind + 3 + layerOffset];
			detectors.zs = z[indz];
			detectors.zd = z[indz + 1];
		}
		else {
			const uint32_t indz = (uint32_t)z_index[oo] * 2U;
			detectors.xs = x[ind];
			detectors.ys = x[ind + 1];
			detectors.xd = x[ind + 2];
			detectors.yd = x[ind + 3];
			detectors.zs = z[indz];
			detectors.zd = z[indz + 1];
		}
	}
	if (n_rays3D > 1)
		multirayCoordinateShiftZ(detectors, lorZ, crz, n_rays3D);
	if (n_rays > 1)
		multirayCoordinateShiftXY(detectors, lorXY, crxy, n_rays);
}

// Transaxial multi-ray shift: the ray is moved perpendicular to the LOR in the xy-plane by ((lor + 0.5) * cr / n_rays - cr / 2),
// both end points by the same amount (GPU: multirayCoordinateShiftXY)
template <typename T>
inline void multirayCoordinateShiftXY(Det<T>& detectors, const int lor, const T cr, const uint16_t n_rays) {
	const T dx = detectors.xd - detectors.xs;
	const T dy = detectors.yd - detectors.ys;
	const T lenXY = std::sqrt(dx * dx + dy * dy);
	if (lenXY < (T)1e-8)
		return;
	const T interval = ((static_cast<T>(lor) + (T)0.5) * (cr / static_cast<T>(n_rays)) - cr * (T)0.5) / lenXY;
	const T shiftX = -dy * interval;
	const T shiftY = dx * interval;
	detectors.xs += shiftX;
	detectors.xd += shiftX;
	detectors.ys += shiftY;
	detectors.yd += shiftY;
}

template <typename T>
inline void multirayCoordinateShiftZ(Det<T>& detectors, const int lor, const T cr, const uint16_t n_rays3D) {
	T interval = cr / (static_cast<T>(n_rays3D * 2));
	detectors.zs += (interval - cr / (T)2.);
	detectors.zd += (interval - cr / (T)2.);
	interval *= (T)2.;
	detectors.zs += interval * lor;
	detectors.zd += interval * lor;
}

// Multi-ray (n_rays > 1 or n_rays3D > 1): the rays are spread uniformly over the detector pixel, i.e. the pixel index is shifted by
// (lor + 0.5) / n_rays - 0.5 (GPU: getDetectorCoordinatesCT). List-mode CT uses the same coordinate shifts as list-mode PET
// (GPU: getDetectorCoordinatesListmode), crxy and crz are then the detector pitches
template <typename T>
inline void get_detector_coordinates_CT(const T* x,const T* z, const uint32_t size_x, Det<T>& detectors, const int64_t lo, const uint32_t subsets,
	const uint32_t size_y, const int64_t ix, const int64_t iy, const int64_t iz, const T dPitch, const int64_t nProjections,
	const uint8_t list_mode_format, const bool pitch, const int lorXY = 0, const int lorZ = 0, const uint16_t n_rays = 1, const uint16_t n_rays3D = 1,
	const T crxy = (T)0., const T crz = (T)0.) {
	if (list_mode_format > 0) {
		detectors.xs = x[lo * 6];
		detectors.ys = x[lo * 6 + 1];
		detectors.zs = x[lo * 6 + 2];
		detectors.xd = x[lo * 6 + 3];
		detectors.yd = x[lo * 6 + 4];
		detectors.zd = x[lo * 6 + 5];
		if (n_rays3D > 1)
			multirayCoordinateShiftZ(detectors, lorZ, crz, n_rays3D);
		if (n_rays > 1)
			multirayCoordinateShiftXY(detectors, lorXY, crxy, n_rays);
	}
	else {
		const int id = iz * 6;
		detectors.xs = x[id];
		detectors.ys = x[id + 1];
		detectors.zs = x[id + 2];
		detectors.xd = x[id + 3];
		detectors.yd = x[id + 4];
		detectors.zd = x[id + 5];
		T indeksi1 = static_cast<T>(ix) - static_cast<T>(size_x) / 2. + .5;
		T indeksi2 = static_cast<T>(iy) - static_cast<T>(size_y) / 2. + .5;
		if (n_rays > 1)
			indeksi1 += (static_cast<T>(lorXY) + (T).5) / static_cast<T>(n_rays) - (T).5;
		if (n_rays3D > 1)
			indeksi2 += (static_cast<T>(lorZ) + (T).5) / static_cast<T>(n_rays3D) - (T).5;
		if (pitch) {
			const int idz = iz * 6;
			detectors.xd += z[idz] * indeksi1 + z[idz + 3] * indeksi2;
			detectors.yd += z[idz + 1] * indeksi1 + z[idz + 4] * indeksi2;
			detectors.zd += z[idz + 2] * indeksi1 + z[idz + 5] * indeksi2;
		}
		else {
			const int idz = iz * 2;
			detectors.xd += z[idz] * indeksi1;
			detectors.yd += z[idz + 1] * indeksi1;
			detectors.zd += dPitch * indeksi2;
		}
	}
}

// Clips the SPECT ray s -> d to the support (GPU: extendRayToEllipse): the end points are moved to where the line through s and d enters
// and exits the support, the direction s -> d is kept. The start point is only moved if it is outside the support. If the line does not
// intersect the support, d is set to s (zero-length ray, skipped by the caller). bmax* are the far edges of the image volume, used as the
// box support if the support has not been set (non-positive radius)
template <typename T>
inline void extendRayToSupport(Det<T>& detectors, const paramStruct<T>& param, const T bmaxx, const T bmaxy, const T bmaxz) {
	const T epsVal = (T)1e-8;
	const T sx = detectors.xs, sy = detectors.ys, sz = detectors.zs;
	const T dirx = detectors.xd - sx, diry = detectors.yd - sy, dirz = detectors.zd - sz;
	const bool supportSet = param.ellipseRadiusX > (T)0. && param.ellipseRadiusY > (T)0. && param.ellipseRadiusZ > (T)0.;
	// Box support (also the fallback when no support has been set). The hosts pass FLT_MAX (or Inf) as the power of the box
	if (!supportSet || param.ellipsePower > (T)1e20) {
		T boxMin[3], boxMax[3];
		if (supportSet) {
			boxMin[0] = param.ellipseCenterX - param.ellipseRadiusX;
			boxMin[1] = param.ellipseCenterY - param.ellipseRadiusY;
			boxMin[2] = param.ellipseCenterZ - param.ellipseRadiusZ;
			boxMax[0] = param.ellipseCenterX + param.ellipseRadiusX;
			boxMax[1] = param.ellipseCenterY + param.ellipseRadiusY;
			boxMax[2] = param.ellipseCenterZ + param.ellipseRadiusZ;
		}
		else {
			boxMin[0] = param.bx;
			boxMin[1] = param.by;
			boxMin[2] = param.bz;
			boxMax[0] = bmaxx;
			boxMax[1] = bmaxy;
			boxMax[2] = bmaxz;
		}
		const T p0[3] = { sx, sy, sz };
		const T dir[3] = { dirx, diry, dirz };
		T tmin = (T)-1e8;
		T tmax = (T)1e8;
		for (int kk = 0; kk < 3; kk++) {
			if (std::fabs(dir[kk]) < epsVal) {
				if (p0[kk] < boxMin[kk] || p0[kk] > boxMax[kk]) {
					detectors.xd = sx;
					detectors.yd = sy;
					detectors.zd = sz;
					return;
				}
			}
			else {
				const T t1 = (boxMin[kk] - p0[kk]) / dir[kk];
				const T t2 = (boxMax[kk] - p0[kk]) / dir[kk];
				tmin = std::max(tmin, std::min(t1, t2));
				tmax = std::min(tmax, std::max(t1, t2));
			}
		}
		if (tmax < tmin) {
			detectors.xd = sx;
			detectors.yd = sy;
			detectors.zd = sz;
			return;
		}
		if (!((sx >= boxMin[0] && sx <= boxMax[0]) && (sy >= boxMin[1] && sy <= boxMax[1]) && (sz >= boxMin[2] && sz <= boxMax[2]))) {
			detectors.xs = sx + tmin * dirx;
			detectors.ys = sy + tmin * diry;
			detectors.zs = sz + tmin * dirz;
		}
		detectors.xd = sx + tmax * dirx;
		detectors.yd = sy + tmax * diry;
		detectors.zd = sz + tmax * dirz;
	}
	// Elliptic cylinder along z
	else if (param.ellipsePower == (T)2.) {
		const T cz = param.ellipseCenterZ;
		const T az = param.ellipseRadiusZ;
		const T px = sx - param.ellipseCenterX;
		const T py = sy - param.ellipseCenterY;
		const T invAx2 = (T)1. / (param.ellipseRadiusX * param.ellipseRadiusX);
		const T invAy2 = (T)1. / (param.ellipseRadiusY * param.ellipseRadiusY);
		const T A = dirx * dirx * invAx2 + diry * diry * invAy2;
		const T B = (T)2. * (px * dirx * invAx2 + py * diry * invAy2);
		const T C = px * px * invAx2 + py * py * invAy2 - (T)1.;
		T tmin = (T)-1e8;
		T tmax = (T)1e8;
		if (A < epsVal) {
			if (C > (T)0.) {
				detectors.xd = sx;
				detectors.yd = sy;
				detectors.zd = sz;
				return;
			}
		}
		else {
			const T discriminant = B * B - (T)4. * A * C;
			if (discriminant < (T)0.) {
				detectors.xd = sx;
				detectors.yd = sy;
				detectors.zd = sz;
				return;
			}
			const T t1 = (-B - std::sqrt(discriminant)) / ((T)2. * A);
			const T t2 = (-B + std::sqrt(discriminant)) / ((T)2. * A);
			tmin = std::max(tmin, t1);
			tmax = std::min(tmax, t2);
		}
		if (std::fabs(dirz) < epsVal) {
			if (sz < cz - az || sz > cz + az) {
				detectors.xd = sx;
				detectors.yd = sy;
				detectors.zd = sz;
				return;
			}
		}
		else {
			const T t1 = (cz - az - sz) / dirz;
			const T t2 = (cz + az - sz) / dirz;
			tmin = std::max(tmin, std::min(t1, t2));
			tmax = std::min(tmax, std::max(t1, t2));
		}
		if (tmax < tmin) {
			detectors.xd = sx;
			detectors.yd = sy;
			detectors.zd = sz;
			return;
		}
		if (C > (T)0. || sz < cz - az || sz > cz + az) {
			detectors.xs = sx + tmin * dirx;
			detectors.ys = sy + tmin * diry;
			detectors.zs = sz + tmin * dirz;
		}
		detectors.xd = sx + tmax * dirx;
		detectors.yd = sy + tmax * diry;
		detectors.zd = sz + tmax * dirz;
	}
	// Other powers: the ray is not clipped (GPU)
}

// SPECT detector coordinates for the ray lor of the detector element (ix, iy) of the projection iz (GPU: getDetectorCoordinatesSPECT).
// x has the start (s) and end (d) points of the central ray of each projection (6 values per projection), z the transaxial direction
// vector of the detector panel (2 values per projection, 6 with pitch). Both end points are moved to the detector element (transaxial
// pitch dPitchXY along the panel, axial pitch dPitchZ along z) and by the collimator ray shifts (mm) of the ray, the element and the
// detector head. The ray shifts have the compact layout [2, nRays, size_x, size_y, nHeads] (x/y-shift along the panel, z-shift). The ray
// is then clipped to the support (extendRayToSupport), the direction s -> d is kept
template <typename T>
inline void get_detector_coordinates_SPECT(const T* x, const T* z, Det<T>& detectors, const paramStruct<T>& param, const uint32_t lor, const int64_t ix,
	const int64_t iy, const int64_t iz, const uint32_t head, const uint32_t nRays, const T bmaxx, const T bmaxy, const T bmaxz) {
	const int64_t idx = iz * 6LL;
	detectors.xs = x[idx + 0];
	detectors.ys = x[idx + 1];
	detectors.zs = x[idx + 2];
	detectors.xd = x[idx + 3];
	detectors.yd = x[idx + 4];
	detectors.zd = x[idx + 5];

	// Shift from the center of the panel to the current detector element
	const T shiftXY = param.dPitchXY * (static_cast<T>(ix) + ((T)1. - static_cast<T>(param.size_x)) * (T)0.5);
	const T shiftZ = param.dPitchZ * (static_cast<T>(iy) + ((T)1. - static_cast<T>(param.size_y)) * (T)0.5);
	// Index of the ray shifts of the current ray, detector element and detector head
	const size_t idShift = 2ULL * static_cast<size_t>(lor) + 2ULL * static_cast<size_t>(nRays) * (static_cast<size_t>(ix) + static_cast<size_t>(iy) * static_cast<size_t>(param.size_x) +
		static_cast<size_t>(head) * static_cast<size_t>(param.size_x) * static_cast<size_t>(param.size_y));
	const T shiftDXY = (param.rayShiftsDetector != nullptr) ? param.rayShiftsDetector[idShift] : (T)0.;
	const T shiftDZ = (param.rayShiftsDetector != nullptr) ? param.rayShiftsDetector[idShift + 1] : (T)0.;
	const T shiftSXY = (param.rayShiftsSource != nullptr) ? param.rayShiftsSource[idShift] : (T)0.;
	const T shiftSZ = (param.rayShiftsSource != nullptr) ? param.rayShiftsSource[idShift + 1] : (T)0.;

	// Direction vector of the detector panel
	const int64_t idz = iz * (param.pitch ? 6LL : 2LL);
	const T apuX = z[idz];
	const T apuY = z[idz + 1];
	detectors.xd += apuX * (shiftXY + shiftDXY);
	detectors.yd += apuY * (shiftXY + shiftDXY);
	detectors.zd += shiftZ + shiftDZ;
	detectors.xs += apuX * (shiftXY + shiftSXY);
	detectors.ys += apuY * (shiftXY + shiftSXY);
	detectors.zs += shiftZ + shiftSZ;

	extendRayToSupport(detectors, param, bmaxx, bmaxy, bmaxz);
}

// Compute the distance that the ray traverses in the current voxel
template <typename T>
inline T voxelValue(const T t, const T tc, const T L) {
	return (t - tc) * L;
}

// Compute the index of the current voxel
inline uint32_t compute_ind(const int tempj, const int tempi, const int tempk, const uint32_t d_Nx, const uint32_t d_Nyx) {
	uint32_t local_ind = (uint32_t)(tempj)*d_Nx + (uint32_t)(tempi)+(uint32_t)(tempk)*d_Nyx;
	return local_ind;
}

// compute the distance that the ray traverses in the current voxel
template <typename T>
inline T compute_element(T& t0, T& tc, const T L, const T tu, const int u, int& temp_ijk) {
	T local_ele = voxelValue(t0, tc, L);
	temp_ijk += u;
	tc = t0;
	t0 += tu;
	return local_ele;
}
void setThreads();

// Exact bin-integrated TOF weights, the same model as in the GPU kernels (general_opencl_functions.h, TOFH ... TOFLoop)
//
// For an emission at the TOF coordinate x (signed distance from the LOR center, measured in the same way as the bin centers
// TOFCenter) the probability of the event being recorded in the bin with the edges [eLo, eHi] is
//     p(x) = Phi((x - eLo) / sigma) - Phi((x - eHi) / sigma),
// where Phi is the standard normal CDF. The weight of a bin is the mean of p(x) over the ray segment (element) that starts at D and
// ends at D + element, i.e. the exact (not sampled) bin probability averaged over the segment. The weights are NOT normalized over
// the bins: the sum over the bins is 1 when the TOF window covers the segment and smaller when it does not.
// With G(u) = u * Phi(u) + phi(u) (G' = Phi) the segment integral is
//     sigma * [G((x1 - eLo) / sigma) - G((x0 - eLo) / sigma) - G((x1 - eHi) / sigma) + G((x0 - eHi) / sigma)],
// G(u) is split into max(u, 0) + H(|u|), H(t) = phi(t) - t * Q(t) (Q = upper tail). The max(u, 0) terms give the clamp terms s0 and s1
// below. H is zero for t >= TOF_CUTOFF. The weights are clamped to [0, 1]. For nBins == 1 all events are accepted (weight 1).
// The bin width is the spacing of the bin centers, |TOFCenter[1] - TOFCenter[0]|.
// The GPU evaluates the NBINS + 1 unique bin edges once each, which requires the center-first bin order (0, -w, +w, -2w, ...), and
// gives NaN otherwise. Here each bin is evaluated from its own center, which is the same value for that order and does not depend
// on the order of the bins.
// The H terms at the start of a segment are the ones at the end of the previous segment, they are cached in TOFF (see TOFInit).

// H(t) = phi(t) - t * Q(t), t >= 0, computed as phi(t) * s^2 * (c0 + c1 s + ... + c6 s^6), s = 1 / (1 + 0.3 t), the same polynomial
// fit as in the GPU kernels (TOFH, maximum absolute error 3.5e-9), exactly zero after TOF_CUTOFF sigmas
template <typename T>
inline T TOFH(const T t) {
	const T s = (T)1. / ((T)1. + (T)0.3 * t);
	T g = s * (T)0.278036759 + (T)-1.48782828;
	g = s * g + (T)2.71388525;
	g = s * g + (T)-1.70602154;
	g = s * g + (T)1.12047934;
	g = s * g + (T)-0.029800149;
	g = s * g + (T)0.111248612;
	return ((t < (T)TOF_CUTOFF) ? (T)TOF_INV_SQRT_2PI * std::exp((T)-0.5 * t * t) : (T)0.) * s * s * g;
}

// The tail part of the antiderivative at the position x for a bin with the edges eLo and eHi
template <typename T>
inline T TOFTail(const T x, const T eLo, const T eHi, const T invSigma) {
	return TOFH(std::fabs((x - eLo) * invSigma)) - TOFH(std::fabs((x - eHi) * invSigma));
}

// Initializes the cache of the H terms at the current position D, has to be called every time D is set (the start of each ray)
template <typename T>
inline void TOFInit(const T D, const T* TOFCenter, const T invSigma, std::vector<T>& TOFF, const uint32_t nBins) {
	if (nBins == 1)
		return;
	const T halfW = (T)0.5 * std::fabs(TOFCenter[1] - TOFCenter[0]);
	for (uint32_t to = 0; to < nBins; to++)
		TOFF[to] = TOFTail(D, TOFCenter[to] - halfW, TOFCenter[to] + halfW, invSigma);
}

// Computes the weights of the bins for the segment starting at D and ending at D + element (the ray is always traversed from the
// source towards the detector, i.e. D increases, GPU: DDsign = -1) and moves D to the end of the segment
template <typename T>
inline void TOFLoop(const T element, const T* TOFCenter, const T invSigma, T& D, std::vector<T>& TOFWeights, std::vector<T>& TOFF, const uint32_t nBins) {
	const T x0 = D;
	const T x1 = x0 + element;
	D = x1;
	if (nBins == 1) {
		TOFWeights[0] = (T)1.;
		return;
	}
	const T halfW = (T)0.5 * std::fabs(TOFCenter[1] - TOFCenter[0]);
	const T cut = (T)TOF_CUTOFF / invSigma;
	const T du = (x1 - x0) * invSigma;
	for (uint32_t to = 0; to < nBins; to++) {
		const T eLo = TOFCenter[to] - halfW;
		const T eHi = TOFCenter[to] + halfW;
		T wt = (T)0.;
		T F1 = (T)0.;
		// Segments farther than TOF_CUTOFF sigmas from the bin have exactly zero H terms and weights
		if (std::max(x0, x1) >= eLo - cut && std::min(x0, x1) <= eHi + cut) {
			F1 = TOFTail(x1, eLo, eHi, invSigma);
			if (du != (T)0.) {
				const T s0 = std::min(std::max(x0, eLo), eHi);
				const T s1 = std::min(std::max(x1, eLo), eHi);
				wt = ((s1 - s0) * invSigma + (F1 - TOFF[to])) / du;
				wt = (wt < (T)0.) ? (T)0. : ((wt > (T)1.) ? (T)1. : wt);
			}
		}
		TOFF[to] = F1;
		TOFWeights[to] = wt;
	}
}

// TOF coordinate (signed distance from the center of the LOR) at the ray parameter tStart (GPU: TOFDis)
template <typename T>
inline void TOFDis(const T tStart, const T L, T& D) {
	D = L * tStart - L / (T)2.;
}

// Compute the sum in the attenuation correction (distance times attenuation coefficient)
template <typename T>
inline void compute_attenuation(const T val, const uint32_t ind, const T* atten, T& jelppi) {
	jelppi += (val * -atten[ind]);
}

// Computes the forward projection
template <typename T>
inline void forwardProject(const T local_ele, T& ax, const uint32_t local_ind, const T* input) {
	ax = (local_ele * input[local_ind]);
}

// Computes the forward projection, TOFWeights are the weights of the current voxel (see TOFLoop)
template <typename T>
inline void denominator(std::vector<T>& ax, const uint32_t local_ind, T local_ele, const T* input, const bool TOF, const std::vector<T>& TOFWeights,
	const uint32_t nBins, const int lor, const uint16_t nRays) {
	T apu = (T)0.;
	forwardProject(local_ele, apu, local_ind, input);
	if (TOF) {
		for (uint32_t to = 0; to < nBins; to++) {
			if (nRays > 1)
				ax[(size_t)to + (size_t)nBins * (size_t)lor] += apu * TOFWeights[to];
			else
				ax[to] += apu * TOFWeights[to];
		}
	}
	else {
		if (nRays > 1)
			ax[lor] += apu;
		else
			ax[0] += apu;
	}
}

// Compute the backprojection, TOFWeights are the weights of the current voxel (see TOFLoop)
template <typename T>
inline void rhs(const T local_ele, const std::vector<T>& ax, const uint32_t local_ind, T* output, const bool no_norm, T* sensImage,
	const std::vector<T>& TOFWeights, const bool TOF, const uint32_t nBins) {
	T yaxTOF = (T)0.;
	T val = (T)0.;
	if (TOF) {
		for (uint32_t to = 0; to < nBins; to++) {
			const T apu = local_ele * TOFWeights[to];
			val += apu;
			yaxTOF += apu * ax[to];
		}
	}
	else {
		yaxTOF = ax[0] * local_ele;
		val = local_ele;
	}
#pragma omp atomic
	output[local_ind] += yaxTOF;
	if (no_norm == false)
#pragma omp atomic
		sensImage[local_ind] += val;
}

template <typename T>
inline void forwardProjectAF(T* output, std::vector<T>& ax, size_t idx, const T temp, const size_t kk, const bool CT) {

	if (!CT)
		output[idx] += ax[kk] * temp;
	else
		output[idx] += ax[kk];
}

// Compute the voxel index where the current perpendicular measurement starts
template <typename T>
inline int perpendicular_start(const T d_b, const T d, const T d_d, const uint32_t d_N) {
	int tempi = 0;
	T start = d_b - d + d_d;
	for (uint32_t ii = 0u; ii < d_N; ii++) {
		if (start > (T)0.) {
			tempi = (int)(ii);
			break;
		}
		start += d_d;
	}
	return tempi;
}

// Voxel index of a coordinate that is constant along the ray, u = (coordinate - b) / voxel size is the coordinate in voxel units:
// floor(u), except that a ray on a boundary between voxels k - 1 and k (closer than 1e-4 voxels to it) is in the voxel k
// (GPU: tieFloor with BOUNDARY_TOL)
template <typename T>
inline int tieFloor(const T u) {
	const T r = std::floor(u + (T)0.5);
	return static_cast<int>((std::fabs(u - r) < (T)1e-4) ? r : std::floor(u));
}

// Compute the voxel index where the current perpendicular measurement starts, zero if the coordinate is outside the image
// (GPU: perpendicular_start)
template <typename T>
inline int perpendicular_start_tie(const T d_b, const T d, const T d_d, const uint32_t d_N) {
	int tempi = tieFloor((d - d_b) / d_d);
	if (tempi < 0 || tempi >= static_cast<int>(d_N))
		tempi = 0;
	return tempi;
}

// Compute the initial and maximum voxel indices and the direction of the ray, detector greater than source
template <typename T>
inline void d_g_s(const T tmin, const T t_min, const T tmax, const T t_max, uint32_t& v_min, uint32_t& v_max, T& t_0, int32_t& v_u,
	const T diff, const T b, const T d, const T s, const uint32_t N) {

	if (tmin == t_min)
		// (11)
		v_min = 1u;
	else {
		// (2) and (19)
		const T p_t = s + tmin * (diff);
		// (12)
		v_min = static_cast<uint32_t>(ceil((p_t - b) / d));
	}
	t_0 += static_cast<T>(v_min) * d / (diff);
	//  (29)
	v_u = 1;
}

// Compute the initial and maximum voxel indices and the direction of the ray, source greater than detector
template <typename T>
inline void s_g_d(const T tmin, const T t_min, const T tmax, const T t_max, uint32_t& v_min, uint32_t& v_max, T& t_0, int32_t& v_u,
	const T diff, const T b, const T d, const T s, const uint32_t N) {

	if (tmin == t_min)
		// (15)
		v_max = N - 1u;
	else {
		// (2) and (19)
		const T p_t = s + tmin * (diff);
		// (16)
		v_max = static_cast<uint32_t>(floor((p_t - b) / d));
	}
	// (9)
	t_0 += static_cast<T>(v_max) * d / (diff);
	// (29)
	v_u = -1;
}

// Compute the initial and maximum voxel indices and the direction of the ray, detector greater than source
template <typename T>
inline void d_g_s_precomp(const T tmin, const T t_min, const T tmax, const T t_max, uint32_t& v_min, uint32_t& v_max, T& t_0,
	int32_t& v_u, const T diff, const T b, const T d, const T s, const uint32_t N) {

	if (tmin == t_min)
		// (11)
		v_min = 1u;
	else {
		// (2) and (19)
		const T p_t = s + tmin * (diff);
		// (12)
		v_min = static_cast<uint32_t>(ceil((p_t - b) / d));
	}
	if (tmax == t_max)
		// (13)
		v_max = N;
	else {
		// (2) and (19)
		const T p_t = s + tmax * (diff);
		// (14)
		v_max = static_cast<uint32_t>(floor((p_t - b) / d));
	}
	// (9)
	t_0 += static_cast<T>(v_min) * d / (diff);
	//  (29)
	v_u = 1;
}

// Compute the initial and maximum voxel indices and the direction of the ray, source greater than detector
template <typename T>
inline void s_g_d_precomp(const T tmin, const T t_min, const T tmax, const double t_max, uint32_t& v_min, uint32_t& v_max, T& t_0,
	int32_t& v_u, const T diff, const T b, const T d, const T s, const uint32_t N) {

	if (tmin == t_min)
		// (15)
		v_max = N - 1u;
	else {
		// (2) and (19)
		const T p_t = s + tmin * (diff);
		// (16)
		v_max = static_cast<uint32_t>(floor((p_t - b) / d));
	}
	if (tmax == t_max)
		// (17)
		v_min = 0u;
	else {
		// (2) and (19)
		const T p_t = s + tmax * (diff);
		// (18)
		v_min = static_cast<uint32_t>(ceil((p_t - b) / d));
	}
	// (9)
	t_0 += static_cast<T>(v_max) * d / (diff);
	// (29)
	v_u = -1;
}

// Compute the first voxel index
template <typename T>
inline int32_t voxel_index(const T pt, const T diff, const T d, const T apu) {
	return static_cast<int32_t>(std::floor((pt * diff - apu) / d));
}

// Compute whether the ray intercepts the FOV (TYPE == 0), compute the voxel indices of the first voxel intercepted, 
// total number of voxels traversed and coefficients needed for Siddon's algorithm (2D)
template <typename T>
inline bool siddon_pre_loop_2D(const T b1, const T b2, const T diff1, const T diff2, const T max1, const T max2,
	const T d1, const T d2, const uint32_t N1, const uint32_t N2, int32_t& temp1, int32_t& temp2, T& t1u, T& t2u, uint32_t& Np,
	const int TYPE, const T ys, const T xs, const T yd, const T xd, T& tc, int32_t& u1, int32_t& u2, T& t10, T& t20, T* tMaxOut = nullptr) {
	// If neither x- nor y-directions are perpendicular
// Correspond to the equations (9) and (10) from reference [1]
	const T apu_tx = b1 - xs;
	const T apu_ty = b2 - ys;
	t10 = (apu_tx) / (diff1);
	t20 = (apu_ty) / (diff2);
	const T txback = (max1 - xs) / (diff1);
	const T tyback = (max2 - ys) / (diff2);

	// Equations (5-8)
	const T txmin = std::min(t10, txback);
	const T txmax = std::max(t10, txback);
	const T tymin = std::min(t20, tyback);
	const T tymax = std::max(t20, tyback);

	// (3-4)
	tc = std::max(txmin, tymin);
	const T tmax = std::min(txmax, tymax);
	// The ray parameter where the ray exits the FOV (required for the in-FOV length of the ray)
	if (tMaxOut != nullptr)
		*tMaxOut = tmax;

	uint32_t imin = 0U, imax = 0U, jmin = 0U, jmax = 0U;

	if (TYPE == 0) {
		// If true, then the ray/LOR does not intersect the pixel space --> continue to the next LOR
		if (tc >= tmax) {
			return true;
		}

		// (11-14)
		if (xs < xd)
			d_g_s_precomp(tc, txmin, tmax, txmax, imin, imax, t10, u1, diff1, b1, d1, xs, N1);
		// (15-18)
		else
			s_g_d_precomp(tc, txmin, tmax, txmax, imin, imax, t10, u1, diff1, b1, d1, xs, N1);

		//Same as above
		if (ys < yd)
			d_g_s_precomp(tc, tymin, tmax, tymax, jmin, jmax, t20, u2, diff2, b2, d2, ys, N2);
		else
			s_g_d_precomp(tc, tymin, tmax, tymax, jmin, jmax, t20, u2, diff2, b2, d2, ys, N2);

		Np = imax + 1u + jmax + 1u - imin - jmin;

		//tc = tmin;
	}
	else {
		// (11-14)
		if (xs < xd)
			d_g_s(tc, txmin, tmax, txmax, imin, imax, t10, u1, diff1, b1, d1, xs, N1);
		// (15-18)
		else
			s_g_d(tc, txmin, tmax, txmax, imin, imax, t10, u1, diff1, b1, d1, xs, N1);
		//Same as above
		if (ys < yd)
			d_g_s(tc, tymin, tmax, tymax, jmin, jmax, t20, u2, diff2, b2, d2, ys, N2);
		else
			s_g_d(tc, tymin, tmax, tymax, jmin, jmax, t20, u2, diff2, b2, d2, ys, N2);
	}

	// (2) and (19)
	const T pt = (((std::min)(t10, t20) + tc) / 2.);

	// (26)
	temp1 = voxel_index(pt, diff1, d1, apu_tx);
	// (27)
	temp2 = voxel_index(pt, diff2, d2, apu_ty);

	if (TYPE == 0) {
		if (temp1 < 0 || temp1 >= N1 || temp2 < 0 || temp2 >= N2)
			return true;
	}

	// (28)
	t1u = d1 / std::fabs(diff1);
	t2u = d2 / std::fabs(diff2);

	return false;
}

// Compute whether the ray intercepts the FOV (TYPE == 0), compute the voxel indices of the first voxel intercepted, 
// total number of voxels traversed and coefficients needed for Siddon's algorithm (3D)
template <typename T>
inline bool siddon_pre_loop_3D(const T bx, const T by, const T bz, const T x_diff, const T y_diff, const T z_diff,
	const T maxxx, const T maxyy, const T bzb, const T dx, const T dy, const T dz,
	const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, int32_t& tempi, int32_t& tempj, int32_t& tempk, T& tyu, T& txu, T& tzu,
	uint32_t& Np, const int TYPE, const Det<T> detectors, T& tc, int32_t& iu, int32_t& ju, int32_t& ku, T& tx0, T& ty0, T& tz0, T* tMaxOut = nullptr) {

	const T apu_tx = bx - detectors.xs;
	const T apu_ty = by - detectors.ys;
	const T apu_tz = bz - detectors.zs;
	tx0 = (apu_tx) / (x_diff);
	ty0 = (apu_ty) / (y_diff);
	tz0 = (apu_tz) / (z_diff);
	const T txback = (maxxx - detectors.xs) / (x_diff);
	const T tyback = (maxyy - detectors.ys) / (y_diff);
	const T tzback = (bzb - detectors.zs) / (z_diff);

	const T txmin = std::min(tx0, txback);
	const T txmax = std::max(tx0, txback);
	const T tymin = std::min(ty0, tyback);
	const T tymax = std::max(ty0, tyback);
	const T tzmin = std::min(tz0, tzback);
	const T tzmax = std::max(tz0, tzback);

	tc = std::max(std::max(txmin, tzmin), tymin);
	const T tmax = std::min(std::min(txmax, tzmax), tymax);
	// The ray parameter where the ray exits the FOV (required for the in-FOV length of the ray)
	if (tMaxOut != nullptr)
		*tMaxOut = tmax;

	uint32_t imin = 0U, imax = 0U, jmin = 0U, jmax = 0U, kmin = 0U, kmax = 0U;

	if (TYPE == 0) {
		if (tc >= tmax) {
			return true;
		}

		if (detectors.xs < detectors.xd)
			d_g_s_precomp(tc, txmin, tmax, txmax, imin, imax, tx0, iu, x_diff, bx, dx, detectors.xs, Nx);
		else
			s_g_d_precomp(tc, txmin, tmax, txmax, imin, imax, tx0, iu, x_diff, bx, dx, detectors.xs, Nx);

		if (detectors.ys < detectors.yd)
			d_g_s_precomp(tc, tymin, tmax, tymax, jmin, jmax, ty0, ju, y_diff, by, dy, detectors.ys, Ny);
		else
			s_g_d_precomp(tc, tymin, tmax, tymax, jmin, jmax, ty0, ju, y_diff, by, dy, detectors.ys, Ny);

		if (detectors.zs < detectors.zd)
			d_g_s_precomp(tc, tzmin, tmax, tzmax, kmin, kmax, tz0, ku, z_diff, bz, dz, detectors.zs, Nz);
		else
			s_g_d_precomp(tc, tzmin, tmax, tzmax, kmin, kmax, tz0, ku, z_diff, bz, dz, detectors.zs, Nz);

		Np = (kmax - kmin + 1) + (jmax - jmin + 1) + (imax - imin + 1);
	}
	else {
		if (detectors.xs < detectors.xd)
			d_g_s(tc, txmin, tmax, txmax, imin, imax, tx0, iu, x_diff, bx, dx, detectors.xs, Nx);
		else
			s_g_d(tc, txmin, tmax, txmax, imin, imax, tx0, iu, x_diff, bx, dx, detectors.xs, Nx);

		if (detectors.ys < detectors.yd)
			d_g_s(tc, tymin, tmax, tymax, jmin, jmax, ty0, ju, y_diff, by, dy, detectors.ys, Ny);
		else
			s_g_d(tc, tymin, tmax, tymax, jmin, jmax, ty0, ju, y_diff, by, dy, detectors.ys, Ny);

		if (detectors.zs < detectors.zd)
			d_g_s(tc, tzmin, tmax, tzmax, kmin, kmax, tz0, ku, z_diff, bz, dz, detectors.zs, Nz);
		else
			s_g_d(tc, tzmin, tmax, tzmax, kmin, kmax, tz0, ku, z_diff, bz, dz, detectors.zs, Nz);
	}

	const T pt = ((std::min(std::min(tz0, ty0), tx0) + tc) / 2.);

	tempi = voxel_index(pt, x_diff, dx, apu_tx);
	tempj = voxel_index(pt, y_diff, dy, apu_ty);
	tempk = voxel_index(pt, z_diff, dz, apu_tz);

	if (TYPE == 0) {
		if (tempi < 0 || tempi >= Nx || tempj < 0 || tempj >= Ny || tempk < 0 || tempk >= Nz)
			return true;
	}

	txu = dx / fabs(x_diff);
	tyu = dy / fabs(y_diff);
	tzu = dz / fabs(z_diff);

	return false;
}

// The GPU (general_opencl_functions.h) versions of the Siddon pre-loop, used by projectorType123Implementation4. Implementation 1
// (improved_siddon_precomputed.cpp, improved_Siddon_algorithm_discard.cpp) keeps the functions above.
//
// The functions above select the start voxel from the floor of the ray midpoint between the entry tc and the first plane crossing
// min(tx0, ty0, tz0). When a plane crossing coincides with the entry point (e.g. a ray entering through the lower z face at a point
// that lies on an x or y voxel boundary), that midpoint collapses onto the entry point and the start voxel is the floor of a value
// that lies exactly on a voxel boundary: rounding then gives either an index outside the image (the whole ray is discarded) or a
// start voxel that is one off from the first crossing time of that axis (the traversal is shifted by one voxel). The GPU instead
// derives the start voxel from the same plane indices as the first crossing times (start_voxel), so they are consistent by
// construction, and clamps the plane indices so that the start voxel is always inside the image.

// Rounding towards plus infinity to an unsigned integer, negative values saturate to zero (GPU: CUINT_rtp)
template <typename T>
inline uint32_t ceilSatU(const T v) {
	const T c = std::ceil(v);
	if (!(c > (T)0.))
		return 0u;
	if (c >= (T)4294967295.)
		return 0xFFFFFFFFu;
	return static_cast<uint32_t>(c);
}

// Rounding towards zero to an unsigned integer, negative values saturate to zero (GPU: CUINT_sat_rtz)
template <typename T>
inline uint32_t truncSatU(const T v) {
	if (!(v > (T)0.))
		return 0u;
	if (v >= (T)4294967295.)
		return 0xFFFFFFFFu;
	return static_cast<uint32_t>(std::floor(v));
}

// Plane indices and the first crossing time, detector greater than source (GPU: d_g_s_precomp). The first plane index is clamped
// to at least 1 and the first crossing time is computed directly from the plane index
template <typename T>
inline void d_g_s_precomp_clamp(const T tmin, const T t_min, const T tmax, const T t_max, uint32_t& v_min, uint32_t& v_max, T& t_0,
	int32_t& v_u, const T diff, const T b, const T d, const T s, const uint32_t N) {
	if (tmin == t_min)
		// (11)
		v_min = 1u;
	else {
		// (2) and (19)
		const T p_t = s + tmin * (diff);
		// (12)
		const uint32_t v = ceilSatU((p_t - b) / d);
		v_min = v < 1u ? 1u : v;
	}
	if (tmax == t_max)
		// (13)
		v_max = N;
	else {
		// (2) and (19)
		const T p_t = s + tmax * (diff);
		// (14)
		v_max = truncSatU((p_t - b) / d);
	}
	// (9)
	t_0 = (b + static_cast<T>(v_min) * d - s) / diff;
	//  (29)
	v_u = 1;
}

// Plane indices and the first crossing time, source greater than detector (GPU: s_g_d_precomp). The first plane index is clamped
// to at most N - 1 and the first crossing time is computed directly from the plane index
template <typename T>
inline void s_g_d_precomp_clamp(const T tmin, const T t_min, const T tmax, const T t_max, uint32_t& v_min, uint32_t& v_max, T& t_0,
	int32_t& v_u, const T diff, const T b, const T d, const T s, const uint32_t N) {
	if (tmin == t_min)
		// (15)
		v_max = N - 1u;
	else {
		// (2) and (19)
		const T p_t = s + tmin * (diff);
		// (16)
		const uint32_t v = truncSatU((p_t - b) / d);
		v_max = v > N - 1u ? N - 1u : v;
	}
	if (tmax == t_max)
		// (17)
		v_min = 0u;
	else {
		// (2) and (19)
		const T p_t = s + tmax * (diff);
		// (18)
		v_min = ceilSatU((p_t - b) / d);
	}
	// (9)
	t_0 = (b + static_cast<T>(v_max) * d - s) / diff;
	// (29)
	v_u = -1;
}

// Initial voxel index (beginning of the ray) from the same plane indices as the first crossing times, so both stay consistent
// under rounding, clamped to the image (GPU: start_voxel)
inline int32_t start_voxel(const int32_t u, const uint32_t v_min, const uint32_t v_max, const uint32_t N) {
	const int32_t v = (u > 0) ? static_cast<int32_t>(v_min) - 1 : static_cast<int32_t>(v_max);
	return v < 0 ? 0 : (v > static_cast<int32_t>(N) - 1 ? static_cast<int32_t>(N) - 1 : v);
}

// The part of the voxel segment [tc, t0] that is inside the line segment t = [0, 1] (zero if the voxel is outside of it) and the move
// to the next voxel (GPU: compute_element without ORTH). A segment that is completely inside one voxel gets its actual length.
// pass is true if and only if the length is positive
template <typename T>
inline T compute_element_clip(T& t0, T& tc, const T L, const T tu, const int u, int& temp_ijk, bool& pass) {
	const T local_ele = std::max(std::min(t0, (T)1.) - std::max(tc, (T)0.), (T)0.) * L;
	pass = local_ele > (T)0.;
	temp_ijk += u;
	tc = t0;
	t0 += tu;
	return local_ele;
}

// Whether the ray intersects the FOV (true if it does NOT), the first voxel, the number of voxels Np and the Siddon coefficients (2D),
// the same as the GPU siddon_pre_loop_2D. tMaxOut is the ray parameter where the ray exits the FOV
template <typename T>
inline bool siddon_pre_loop_2D_clamp(const T b1, const T b2, const T diff1, const T diff2, const T max1, const T max2,
	const T d1, const T d2, const uint32_t N1, const uint32_t N2, int32_t& temp1, int32_t& temp2, T& t1u, T& t2u, uint32_t& Np,
	const T ys, const T xs, const T yd, const T xd, T& tc, int32_t& u1, int32_t& u2, T& t10, T& t20, T& tMaxOut) {
	// Correspond to the equations (9) and (10) from reference [1]
	const T apu_tx = b1 - xs;
	const T apu_ty = b2 - ys;
	t10 = (apu_tx) / (diff1);
	t20 = (apu_ty) / (diff2);
	const T txback = (max1 - xs) / (diff1);
	const T tyback = (max2 - ys) / (diff2);

	// Equations (5-8)
	const T txmin = std::min(t10, txback);
	const T txmax = std::max(t10, txback);
	const T tymin = std::min(t20, tyback);
	const T tymax = std::max(t20, tyback);

	// (3-4)
	tc = std::max(txmin, tymin);
	const T tmax = std::min(txmax, tymax);
	tMaxOut = tmax;

	// If true, then the ray/LOR does not intersect the pixel space --> continue to the next LOR
	if (tc >= tmax)
		return true;

	uint32_t imin = 0U, imax = 0U, jmin = 0U, jmax = 0U;
	// (11-14)
	if (xs < xd)
		d_g_s_precomp_clamp(tc, txmin, tmax, txmax, imin, imax, t10, u1, diff1, b1, d1, xs, N1);
	// (15-18)
	else
		s_g_d_precomp_clamp(tc, txmin, tmax, txmax, imin, imax, t10, u1, diff1, b1, d1, xs, N1);
	//Same as above
	if (ys < yd)
		d_g_s_precomp_clamp(tc, tymin, tmax, tymax, jmin, jmax, t20, u2, diff2, b2, d2, ys, N2);
	else
		s_g_d_precomp_clamp(tc, tymin, tmax, tymax, jmin, jmax, t20, u2, diff2, b2, d2, ys, N2);

	Np = imax + 1u + jmax + 1u - imin - jmin;

	// (26)
	temp1 = start_voxel(u1, imin, imax, N1);
	// (27)
	temp2 = start_voxel(u2, jmin, jmax, N2);

	// (28)
	t1u = d1 / std::fabs(diff1);
	t2u = d2 / std::fabs(diff2);

	return false;
}

// Same as above, but for the 3D case (GPU: siddon_pre_loop_3D)
template <typename T>
inline bool siddon_pre_loop_3D_clamp(const T bx, const T by, const T bz, const T x_diff, const T y_diff, const T z_diff,
	const T maxxx, const T maxyy, const T bzb, const T dx, const T dy, const T dz,
	const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, int32_t& tempi, int32_t& tempj, int32_t& tempk, T& txu, T& tyu, T& tzu,
	uint32_t& Np, const Det<T>& detectors, T& tc, int32_t& iu, int32_t& ju, int32_t& ku, T& tx0, T& ty0, T& tz0, T& tMaxOut) {

	tx0 = (bx - detectors.xs) / (x_diff);
	ty0 = (by - detectors.ys) / (y_diff);
	tz0 = (bz - detectors.zs) / (z_diff);
	const T txback = (maxxx - detectors.xs) / (x_diff);
	const T tyback = (maxyy - detectors.ys) / (y_diff);
	const T tzback = (bzb - detectors.zs) / (z_diff);

	const T txmin = std::min(tx0, txback);
	const T txmax = std::max(tx0, txback);
	const T tymin = std::min(ty0, tyback);
	const T tymax = std::max(ty0, tyback);
	const T tzmin = std::min(tz0, tzback);
	const T tzmax = std::max(tz0, tzback);

	tc = std::max(std::max(txmin, tzmin), tymin);
	const T tmax = std::min(std::min(txmax, tzmax), tymax);
	tMaxOut = tmax;

	if (tc >= tmax)
		return true;

	uint32_t imin = 0U, imax = 0U, jmin = 0U, jmax = 0U, kmin = 0U, kmax = 0U;

	if (detectors.xs < detectors.xd)
		d_g_s_precomp_clamp(tc, txmin, tmax, txmax, imin, imax, tx0, iu, x_diff, bx, dx, detectors.xs, Nx);
	else
		s_g_d_precomp_clamp(tc, txmin, tmax, txmax, imin, imax, tx0, iu, x_diff, bx, dx, detectors.xs, Nx);

	if (detectors.ys < detectors.yd)
		d_g_s_precomp_clamp(tc, tymin, tmax, tymax, jmin, jmax, ty0, ju, y_diff, by, dy, detectors.ys, Ny);
	else
		s_g_d_precomp_clamp(tc, tymin, tmax, tymax, jmin, jmax, ty0, ju, y_diff, by, dy, detectors.ys, Ny);

	if (detectors.zs < detectors.zd)
		d_g_s_precomp_clamp(tc, tzmin, tmax, tzmax, kmin, kmax, tz0, ku, z_diff, bz, dz, detectors.zs, Nz);
	else
		s_g_d_precomp_clamp(tc, tzmin, tmax, tzmax, kmin, kmax, tz0, ku, z_diff, bz, dz, detectors.zs, Nz);

	Np = (kmax - kmin + 1u) + (jmax - jmin + 1u) + (imax - imin + 1u);

	tempi = start_voxel(iu, imin, imax, Nx);
	tempj = start_voxel(ju, jmin, jmax, Ny);
	tempk = start_voxel(ku, kmin, kmax, Nz);

	txu = dx / std::fabs(x_diff);
	tyu = dy / std::fabs(y_diff);
	tzu = dz / std::fabs(z_diff);

	return false;
}

// Compute the probability for one emission in perpendicular detector case
template <typename T>
inline T perpendicular_elements(const uint32_t N, const T dd, const T dd2, const T b, const T d, const uint32_t N1, const uint32_t N2, const T* atten,
	const T local_norm, const bool attenuation_correction, const bool normalization, const bool CTAttenuation, 
	int32_t& tempk, const uint32_t NN, const T global_factor, const bool scatter, const T local_scat, int& indX, int& indY,
	const int indZ, const T L, const uint16_t nRays, const int projType, const uint32_t Nx, const uint32_t Ny, const bool CT, const int64_t idx, 
	const bool SPECT = false) {
	uint32_t apu = 0u;
	// Find the closest y-index value by finding the smallest y-distance between detector 2 and all the y-pixel coordinates
	T start = b - dd + d;
	for (int32_t ii = 0u; ii < N1; ii++) {
		if (start > 0.) {
			apu = (ii);
			break;
		}
		start += d;
	}
	tempk = apu * N + tempk * N1 * N2;
	if (N == 1)
		indX = apu;
	else
		indY = apu;
	if (!CT) {
		T temp = 1.;
		if (nRays > 1)
			temp = (T)1. / (L * static_cast<T>(nRays));
		else if (projType == 1)
			temp = (T)1. / L;
		if (attenuation_correction && CTAttenuation && !SPECT) {
			T jelppi = 0.;
			uint32_t atnindX = indX;
			uint32_t atnindY = indY;
			uint32_t atnindZ = indZ;
			for (uint32_t iii = 0u; iii < N2; iii++) {
				if (NN == 1)
					atnindX = iii;
				else
					atnindY = iii;
				jelppi += (dd2 * -atten[atnindX + atnindY * Nx + atnindZ * Nx * Ny]);
			}
			temp *= std::exp(jelppi);
		}
		if (normalization)
			temp *= local_norm;
		if (!CTAttenuation && attenuation_correction)
			temp *= atten[idx];
		if (scatter)
			temp *= local_scat;
		temp *= global_factor;

		return temp;
	}
	else
		if (nRays > 1)
			return 1. / static_cast<T>(nRays);
		else
			return 1.;
}

// Compute the probability for one emission in the axis-parallel case (projectorType123Implementation4), the same as the GPU
// perpendicular_elements:
// N is the stride of the perpendicular coordinate, NN the stride of the coordinate along the ray, N1 the number of voxels in the
// perpendicular direction and N1 * N2 = Nx * Ny. b, d and dd are the edge of the image, the voxel size and the detector
// coordinate in the perpendicular direction, dd2 the voxel size along the ray. perpY is true if the perpendicular coordinate is y.
// On input tempk is the slice index, on output the linear index of the first voxel (coordinate along the ray at zero).
// The voxels apuX1 ... apuX2 are traversed, dT1 and dT2 are the intersection lengths of the end voxels (see the caller) and
// pathLength the length of the ray inside the FOV
template <typename T>
inline T perpendicular_elements_par(const uint32_t N, const T dd, const T dd2, const T b, const T d, const uint32_t N1, const uint32_t N2, const T* atten,
	const T local_norm, const bool attenuation_correction, const bool normalization, const bool CTAttenuation, int32_t& tempk, const uint32_t NN,
	const T global_factor, const bool scatter, const T local_scat, int& indX, int& indY, const T L, const uint16_t nRays, const int projType,
	const bool CT, const int64_t idx, const bool SPECT, const bool perpY, const T pathLength, const bool useTotLength, const int32_t apuX1,
	const int32_t apuX2, const T dT1, const T dT2) {
	const int32_t apu = perpendicular_start_tie(b, dd, d, N1);
	tempk = apu * N + tempk * N1 * N2;
	if (perpY)
		indY = apu;
	else
		indX = apu;
	if (!CT) {
		// The total length of the ray (useTotLength, GPU -DTOTLENGTH) or the length inside the FOV. For SPECT the total length is the
		// length inside the support (GPU)
		const T LN = useTotLength ? L : pathLength;
		T temp = 1.;
		if (nRays > 1)
			temp = (T)1. / (LN * static_cast<T>(nRays));
		else if (projType == 1)
			temp = (T)1. / LN;
		if (attenuation_correction && CTAttenuation && !SPECT) {
			T jelppi = 0.;
			for (int32_t iii = apuX1; iii <= apuX2; iii++) {
				T d_in = dd2;
				if (iii == apuX1 && iii == apuX2)
					d_in = std::max(dT1 + dT2 - dd2, (T)0.);
				else if (iii == apuX1)
					d_in = dT1;
				else if (iii == apuX2)
					d_in = dT2;
				jelppi += (d_in * -atten[static_cast<uint32_t>(tempk) + static_cast<uint32_t>(iii) * NN]);
			}
			temp *= std::exp(jelppi);
		}
		if (normalization)
			temp *= local_norm;
		if (scatter)
			temp *= local_scat;
		// Measurement-based attenuation (GPU -DATNM), atten is then the attenuation of the measurements of the current subset
		if (!CTAttenuation && attenuation_correction)
			temp *= atten[idx];
		temp *= global_factor;

		return temp;
	}
	else
		if (nRays > 1)
			return 1. / static_cast<T>(nRays);
		else
			return 1.;
}

// this function was taken from: https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
// Currently unused
template<typename T>
inline std::vector<size_t> sort_indexes(const std::vector<T>& v)
{
	// initialize original index locations
	std::vector<size_t> idx(v.size());
	std::iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

	return idx;
}

// Source: http://www.alecjacobson.com/weblog/?p=4544 & https://ideone.com/Z7zldb
// Unused, replaced with OpenMP
//class ThreadPool {
//
//public:
//
//	template<typename Index, typename Callable>
//	static void ParallelFor(Index start, Index end, Callable func) {
//		// Estimate number of threads in the pool
//		const static unsigned nb_threads_hint = std::thread::hardware_concurrency();
//		const static unsigned nb_threads = (nb_threads_hint == 0u ? 8u : nb_threads_hint);
//
//		// Size of a slice for the range functions
//		Index n = end - start + 1;
//		Index slice = (Index)std::round(n / static_cast<double> (nb_threads));
//		slice = std::max(slice, Index(1));
//
//		// [Helper] Inner loop
//		auto launchRange = [&func](uint32_t k1, uint32_t k2) {
//			for (Index k = k1; k < k2; k++) {
//				func(k);
//			}
//		};
//
//		// Create pool and launch jobs
//		std::vector<std::thread> pool;
//		pool.reserve(nb_threads);
//		Index i1 = start;
//		Index i2 = std::min(start + slice, end);
//		for (unsigned i = 0; i + 1 < nb_threads && i1 < end; ++i) {
//			pool.emplace_back(launchRange, i1, i2);
//			i1 = i2;
//			i2 = std::min(i2 + slice, end);
//		}
//		if (i1 < end) {
//			pool.emplace_back(launchRange, i1, end);
//		}
//
//		// Wait for jobs to finish
//		for (std::thread &t : pool) {
//			if (t.joinable()) {
//				t.join();
//			}
//		}
//	}
//
//	// Serial version for easy comparison
//	template<typename Index, typename Callable>
//	static void SequentialFor(Index start, Index end, Callable func) {
//		for (Index i = start; i < end; i++) {
//			func(i);
//		}
//	}
//
//};

//#undef _OPENMP

/// <summary>
/// Compute the forward or backward projection for projector types 1, 2 or 3. The input/output data can be either float or double.
/// </summary>
/// <param name="param most required and optional input parameters. Modify the struct accordingly"></param>
/// <param name="nMeas Number of elements in the input data, i.e. the number of measurements in the current subset"></param>
/// <param name="output pointer to the double or float output array (either the result of the forward projection or backprojection, make sure the dimensions are correct)"></param>
/// <param name="x Transaxial detector coordinates (double or float). The order is x-source, y-source, x-detector, y-detector, i.e. the dimensions should be 4 * total detector elements"></param>
/// <param name="z Axial detector coordinates (double or float). The order is z-source, z-detector, i.e. the dimensions should be 2 * total detector elements"></param>
/// <param name="input pointer to the double or float input array (either the input to the forward projection or backprojection, make sure the dimensions are correct)"></param>
/// <param name="fp if fp == 1, forward projection is computed, if fp == 2 bacprojection. Default is forward projection (optional)"></param>
/// <param name="SensImage (double or float) pointer for the sensitivity image (optional)"></param>
/// <param name="detIndex detector indices for each measurement index for the raw data only. Input this only if you use raw data!"></param>
/// <param name="nCores Number of cores/threads used. You can optionally input the number of threads/cores you want to use. Default uses all threads."></param>
/// <returns></returns>
template <typename T>
void projectorType123Implementation4(paramStruct<T>& param, const int64_t nMeas, T* output, const T* x, const T* z,
	T* input, const bool CT = false, const bool SPECT = false, const uint8_t fp = 1, T* SensImage = nullptr, const uint16_t* detIndex = nullptr, 
	const uint32_t nCores = 0) {

#ifdef _OPENMP
	if (nCores == 0U)
		setThreads();
	else
		omp_set_num_threads(nCores);
#endif

	// Size of a single 2D image
	const uint32_t Nyx = param.Ny * param.Nx;

	//uint32_t bin = 0U;

	// Distance of the last slice from origin
	const T bmaxz = param.bz + static_cast<T>(param.Nz) * param.dz;
	const T bmaxy = static_cast<T>(param.Ny) * param.dy + param.by;
	const T bmaxx = static_cast<T>(param.Nx) * param.dx + param.bx;

	int64_t lo = 0LL;

	uint32_t nRays = param.nRays2D * param.nRays3D;
	const T invSigma = (T)1. / param.sigma_x;
#ifdef _OPENMP
#pragma omp parallel
	{
		std::vector<T> ax(param.nBins);
		// TOF weights of the current voxel and the cached H terms (see TOFLoop)
		std::vector<T> TOFWeights(param.nBins), TOFF(param.nBins);
#if _OPENMP >= 201511 && defined(MATLAB)
#pragma omp for schedule(monotonic:dynamic, nChunks)
#else
#pragma omp for schedule(dynamic, nChunks)
#endif
#else
	std::vector<T> ax(param.nBins);
	std::vector<T> TOFWeights(param.nBins), TOFF(param.nBins);
#endif
	for (int64_t lo = 0LL; lo < nMeas; lo++) {

		int64_t ix = lo, iy = 0, iz = 0, izD = 0;

		if (param.subsets > 1 && param.subsetType == 1 && param.listMode == 0) {
			ix *= param.subsets;
			ix += param.currentSubset;
			iy = ix % param.size_y;
			iz = ix / ((int64_t)param.size_x * (int64_t)param.size_y);
			ix /= param.size_y;
			ix = ix % param.size_x;
		}
		else if (param.subsets > 1 && param.subsetType == 2 && param.listMode == 0) {
			ix *= param.subsets;
			ix += param.currentSubset;
			iz = ix / ((int64_t)param.size_x * (int64_t)param.size_y);
			iy = ix / param.size_x;
			iy = iy % param.size_y;
			ix = ix % param.size_x;
		}
		else if (param.subsets > 1 && param.subsetType == 4 && param.listMode == 0) {
			ix = (ix / param.size_x) * (int64_t)param.size_x * (int64_t)param.subsets + (int64_t)param.size_x * (int64_t)param.currentSubset + ix % param.size_x;
			iz = ix / ((int64_t)param.size_x * (int64_t)param.size_y);
			iy = (ix / param.size_x) % param.size_y;
			ix = ix % param.size_x;
		}
		else if (param.subsets > 1 && param.subsetType == 5 && param.listMode == 0) {
			iy = ix % param.size_y;
			ix = (ix / param.size_y * param.subsets + param.currentSubset);
			iz = ix / param.size_x;
			ix = ix % param.size_x;
		}
		else if ((param.subsetType >= 8 || param.subsets == 1) && param.listMode == 0) {
			iz = lo / ((int64_t)param.size_x * (int64_t)param.size_y);
			ix = (lo - iz * (int64_t)param.size_x * (int64_t)param.size_y) % param.size_x;
			iy = (lo - iz * (int64_t)param.size_x * (int64_t)param.size_y) / param.size_x;
			if (param.subsetType >= 8 && param.subsets > 1) {
				iz += (int64_t)param.nMeas;
			}
		}
		else {
			if ((param.listMode > 0 && param.computeSensIm == 0) || (!(CT || SPECT) && param.listMode == 0)) {
				ix = lo;
				iy = 0;
				iz = 0;
			}
			else if (param.listMode > 0 && param.computeSensIm) {
				// The same enumeration as the GPU (projectorType123.cl, SENS): the pair index is
				// lo = a + b * det_per_ring + (ringA * rings + ringB) * det_per_ring^2 + (layerA * nLayers + layerB) * det_per_ring^2 * rings^2,
				// i.e. nLayers^2 * det_per_ring^2 * rings^2 ordered pairs. The first detector (a, ringA) is the source (zs = z[ringA]), the
				// second (b, ringB) the detector. The detectors of layer l are at the offset l * det_per_ring and the rings at l * rings in the
				// coordinate arrays (the normalization/scatter index is lo itself, as idx on the GPU)
				const int64_t detPerRing = static_cast<int64_t>(param.det_per_ring);
				const int64_t nRings = static_cast<int64_t>(param.rings);
				const int64_t nLayersS = static_cast<int64_t>(param.nLayerCount);
				const int64_t ringPairs = nRings * nRings;
				ix = lo % detPerRing;
				iy = (lo / detPerRing) % detPerRing;
				int64_t ringPair = lo / (detPerRing * detPerRing);
				const int64_t layerPair = ringPair / ringPairs;
				const int64_t layerA = layerPair / nLayersS;
				const int64_t layerB = layerPair % nLayersS;
				ringPair -= layerPair * ringPairs;
				iz = ringPair / nRings + layerA * nRings;
				izD = ringPair % nRings + layerB * nRings;
				ix += layerA * detPerRing;
				iy += layerB * detPerRing;
				// The enumeration contains ordered detector pairs, i.e. every LOR {(detector ix, ring iz), (detector iy, ring izD)} twice, as (A,B) and (B,A).
				// Only one orientation is kept so that every LOR contributes exactly once to the sensitivity image (ix, iy, iz and izD are
				// layer-global indices, ix == iy can only occur within the same layer)
				if (layerPair >= nLayersS * nLayersS || ix > iy || (ix == iy && iz >= izD))
					continue;
			}
			else {
				iz = ix / ((int64_t)param.size_x * (int64_t)param.size_y);
				iy = ix % param.size_y;
				ix = (ix / param.size_y) % param.size_x;
			}
		}
		int64_t ind = lo;
		// Detector head of the current SPECT projection (GPU: d_detectorVector[i.z])
		uint32_t head = 0u;
		if (SPECT && param.detectorVector != nullptr)
			head = param.detectorVector[iz];
		if (param.useMaskFP) {
			const int64_t idx1 = ix + iy * (int64_t)param.size_x;
			const bool maskVal = param.maskFP[idx1];
			if (maskVal == false)
				continue;
		}
		Det<T> detectors; 

		//std::vector<T> ax(param.nBins);
		if (fp == 1)
			ax.resize((size_t)param.nBins * (size_t)nRays);
		std::fill(ax.begin(), ax.end(), (T)0.);

		if (fp == 2 && param.listMode <= 1 && !param.computeSensIm) {
			for (int64_t to = 0LL; to < param.nBins; to++)
				ax[to] = input[lo + to * nMeas];
		}
		else if (param.computeSensIm && param.listMode > 0) {
			for (int64_t to = 0LL; to < param.nBins; to++)
				ax[to] = (T)1.;
		}
		else if (fp == 2) {
			for (int64_t to = 0LL; to < param.nBins; to++)
				ax[to] = input[lo + to * nMeas];
		}
		T local_norm = (T)0.;
		T local_scat = (T)0.;
		if (param.normalizationCorrection)
			local_norm = static_cast<T>(param.normCoefs[lo]);
		if (param.scatterCorrectionMult)
			local_scat = static_cast<T>(param.scatterCoefs[lo]);
		int lor = -1;

		// Loop through the rays
		for (int lorZ = 0; lorZ < param.nRays3D; lorZ++) {
			for (int lorXY = 0; lorXY < param.nRays2D; lorXY++) {
				lor++;

				//if (!SPECT) {
				//Det<T> detectors;
				detectors.xd = (T)0., detectors.xs = (T)0., detectors.yd = (T)0., detectors.ys = (T)0., detectors.zd = (T)0., detectors.zs = (T)0.;
				if (CT) { // CT data
					get_detector_coordinates_CT(x, z, param.size_x, detectors, lo, param.subsets, param.size_y, ix, iy, iz, param.dPitchZ, param.nProjections, param.listMode, param.pitch,
						lorXY, lorZ, param.nRays2D, param.nRays3D, param.dPitchXY, param.dPitchZ);
				} else if (SPECT) { // SPECT data
					get_detector_coordinates_SPECT(x, z, detectors, param, static_cast<uint32_t>(lor), ix, iy, iz, head, nRays, bmaxx, bmaxy, bmaxz);
				} else {
					// Raw data
					// Pure list-mode format (e.g. event-by-event)
					if (param.raw || param.listMode > 0) {
						get_detector_coordinates_raw(x, z, detectors, lo, param.listMode, param.nRays2D, param.nRays3D, ix, iy, iz, izD, lorZ, lorXY, param.dPitchZ, param.dPitchXY, param.computeSensIm);
					}
					// Sinogram data
					else {
						get_detector_coordinates(x, z, param.size_x, param.size_y, detectors, param.xy_index, param.z_index, lo, param.subsetType, param.subsets, ix, iy, iz, param.nRays2D, param.nRays3D, lorZ,
							lorXY, param.dPitchZ, param.dPitchXY, param.nLayers);
					}
				}

				// Calculate the x, y and z distances of the detector pair
				T y_diff = (detectors.yd - detectors.ys);
				T x_diff = (detectors.xd - detectors.xs);
				T z_diff = (detectors.zd - detectors.zs);
				// Skip certain cases (e.g. if the x- and y-coordinates are the same for both detectors, LOR between detector n and n)
				if ((y_diff == 0. && x_diff == 0.) || std::isinf(y_diff) || std::isinf(x_diff) || std::isinf(z_diff) || std::isnan(y_diff) || std::isnan(x_diff) || std::isnan(z_diff)) {
					continue;
				}

				// Number of voxels the ray traverses
				uint32_t Np = 0U;

				T jelppi = (T)0.;
				T temp = (T)1.;
				uint32_t d_N0 = param.Nx;
				uint32_t d_N1 = param.Ny;
				uint32_t d_N2 = 1u;
				uint32_t d_N3 = param.Nx;

				int tempi = 0, tempj = 0, tempk = 0, ux = 0, uy = 0, uz = 0;

				T L = norm(x_diff, y_diff, z_diff);
				// Near-parallel-to-an-axis threshold, relative to the ray length (GPU: parTol). Near-axis rays use the axis-parallel branch
				const T parTol = std::max((T)1e-6, (T)4e-7 * L);
				uint32_t local_ind = 0u;
				int localIndX = 0, localIndY = 0, localIndZ = 0;
				T D = 0.;
				T local_ele = 0.;
				uint8_t maskVal = 1;

				if (std::fabs(z_diff) < (T)1e-6 && (std::fabs(y_diff) < parTol || std::fabs(x_diff) < parTol)) {

					// Ring number
					int32_t tempk = tieFloor((detectors.zs - param.bz) / param.dz);
					if (tempk < 0 || tempk >= static_cast<int32_t>(param.Nz))
						continue;
					int indO = 0;
					T d_b, dd, d_db, d_d2;
					// The coordinate along the ray: difference, edge of the image and source coordinate
					T diff1 = (T)0., b1 = (T)0., s1 = (T)0.;
					bool perpY = true;

					// Detectors are perpendicular
					// Siddon cannot be applied --> trivial to compute
					// apuX1 and apuX2 are the first and last voxels (along the ray) inside the segment between the end points,
					// dT1 and dT2 the lengths of the segment inside them
					int32_t apuX1, apuX2;
					T dT1, dT2;
					if (std::fabs(y_diff) < parTol && detectors.yd <= bmaxy && detectors.yd >= param.by && detectors.ys <= bmaxy && detectors.ys >= param.by) {
						apuX1 = 0;
						apuX2 = static_cast<int32_t>(param.Nx) - 1;
						dT1 = param.dx;
						dT2 = param.dx;
						diff1 = x_diff;
						b1 = param.bx;
						s1 = detectors.xs;
						T dist1, dist2 = (T)0.;
						if (detectors.xs > detectors.xd) {
							dist1 = (param.bx - detectors.xd);
							dist2 = (param.bx + static_cast<T>(param.Nx) * param.dx - detectors.xs);
						}
						else {
							dist1 = (param.bx - detectors.xs);
							dist2 = (param.bx + static_cast<T>(param.Nx) * param.dx - detectors.xd);
						}
						// Both endpoints on the same side of the image: the segment does not intersect it
						if (std::max(detectors.xs, detectors.xd) <= param.bx || std::min(detectors.xs, detectors.xd) >= param.bx + static_cast<T>(param.Nx) * param.dx)
							continue;
						if (dist1 < (T)0.) {
							const int kk = static_cast<int>(std::ceil(-dist1 / param.dx)) - 1;
							if (kk < static_cast<int>(param.Nx)) {
								apuX1 = kk;
								dT1 = std::min(dist1 + static_cast<T>(kk + 1) * param.dx, param.dx);
							}
						}
						if (dist2 > (T)0.) {
							const int mm = static_cast<int>(std::ceil(dist2 / param.dx));
							const int kk = static_cast<int>(param.Nx) - mm;
							if (kk >= apuX1) {
								apuX2 = kk;
								dT2 = std::min(static_cast<T>(mm) * param.dx - dist2, param.dx);
							}
						}
						d_b = param.by;
						dd = detectors.yd;
						d_db = param.dy;
						d_d2 = param.dx;
						d_N0 = param.Ny;
						d_N1 = param.Nx;
						// d_N2/d_N3 are the strides of the perpendicular (y) and the marching (x) coordinate, i.e. the y stride is Nx
						d_N2 = param.Nx;
						d_N3 = 1u;
						perpY = true;
					}
					else if (std::fabs(x_diff) < parTol && detectors.xd <= bmaxx && detectors.xd >= param.bx && detectors.xs <= bmaxx && detectors.xs >= param.bx) {
						apuX1 = 0;
						apuX2 = static_cast<int32_t>(param.Ny) - 1;
						dT1 = param.dy;
						dT2 = param.dy;
						diff1 = y_diff;
						b1 = param.by;
						s1 = detectors.ys;
						T dist1, dist2 = (T)0.;
						if (detectors.ys > detectors.yd) {
							dist1 = (param.by - detectors.yd);
							dist2 = (param.by + static_cast<T>(param.Ny) * param.dy - detectors.ys);
						}
						else {
							dist1 = (param.by - detectors.ys);
							dist2 = (param.by + static_cast<T>(param.Ny) * param.dy - detectors.yd);
						}
						// Both endpoints on the same side of the image: the segment does not intersect it
						if (std::max(detectors.ys, detectors.yd) <= param.by || std::min(detectors.ys, detectors.yd) >= param.by + static_cast<T>(param.Ny) * param.dy)
							continue;
						if (dist1 < (T)0.) {
							const int kk = static_cast<int>(std::ceil(-dist1 / param.dy)) - 1;
							if (kk < static_cast<int>(param.Ny)) {
								apuX1 = kk;
								dT1 = std::min(dist1 + static_cast<T>(kk + 1) * param.dy, param.dy);
							}
						}
						if (dist2 > (T)0.) {
							const int mm = static_cast<int>(std::ceil(dist2 / param.dy));
							const int kk = static_cast<int>(param.Ny) - mm;
							if (kk >= apuX1) {
								apuX2 = kk;
								dT2 = std::min(static_cast<T>(mm) * param.dy - dist2, param.dy);
							}
						}
						d_b = param.bx;
						dd = detectors.xd;
						d_db = param.dx;
						d_d2 = param.dy;
						// Default strides: x is the perpendicular coordinate (stride 1), y the marching coordinate (stride Nx)
						perpY = false;
					}
					else
						continue;
					// The length of the ray inside the FOV
					T rayLength = (apuX1 == apuX2) ? (dT1 + dT2 - d_d2) : (dT1 + dT2 + static_cast<T>(apuX2 - apuX1 - 1) * d_d2);
					if (rayLength <= (T)0.)
						rayLength = d_d2;
					localIndZ = tempk;
					temp = perpendicular_elements_par(d_N2, dd, d_d2, d_b, d_db, d_N0, d_N1, param.atten, local_norm, param.attenuationCorrection, param.normalizationCorrection,
						param.CTAttenuation, tempk, d_N3, param.globalFactor, param.scatterCorrectionMult, local_scat, localIndX, localIndY, L, nRays,
						param.projType, CT, lo, SPECT, perpY, rayLength, param.useTotLength, apuX1, apuX2, dT1, dT2);
					// Linear index of the first voxel of the row, the coordinate along the ray is added below
					const uint32_t rowInd = static_cast<uint32_t>(tempk);
					if (fp == 2 && param.useMaskBP) {
						if (perpY)
							indO = localIndY;
						else
							indO = localIndX;
					}
					// The loop is traversed in the direction of the ray (GPU), so that the TOF coordinate increases from the source and the
					// SPECT attenuation accumulates from the start point of the ray
					const int stepDir = (diff1 < (T)0.) ? -1 : 1;
					const int iiStart = (stepDir > 0) ? apuX1 : apuX2;
					const int iiEnd = (stepDir > 0) ? apuX2 : apuX1;
					if (param.TOF) {
						// The TOF coordinate at the entry of the first voxel, D always increases
						// In PET cases the source and/or detector should never be inside the FOV, but this is still taken into account here
						const T entry = b1 + static_cast<T>((stepDir > 0) ? iiStart : iiStart + 1) * d_d2;
						D = std::max((entry - s1) * static_cast<T>(stepDir), (T)0.) - L * (T)0.5;
						TOFInit(D, param.TOFCenters, invSigma, TOFF, param.nBins);
					}
					const T dFirst = (iiStart == iiEnd) ? std::max(dT1 + dT2 - d_d2, (T)0.) : ((stepDir > 0) ? dT1 : dT2);
					const T dLast = (stepDir > 0) ? dT2 : dT1;

					for (int ii = iiStart; ; ii += stepDir) {
						T d_in = d_d2;
						if (ii == iiStart)
							d_in = dFirst;
						else if (ii == iiEnd)
							d_in = dLast;
						const T dTOF = d_in;
						local_ind = rowInd + static_cast<uint32_t>(ii) * d_N3;
						if (param.attenuationCorrection && SPECT && param.CTAttenuation) {
							compute_attenuation(d_in, local_ind, param.atten, jelppi);
							if (param.projType == 1)
								d_in *= std::exp(jelppi);
						}
						if (param.TOF)
							TOFLoop(dTOF, param.TOFCenters, invSigma, D, TOFWeights, TOFF, param.nBins);
						if (fp == 1) {
							denominator(ax, local_ind, d_in, input, param.TOF, TOFWeights, param.nBins, lor, nRays);
						}
						else if (fp == 2) {
							if (param.useMaskBP) {
								maskVal = param.maskBP[indO * d_N2 + ii * d_N3];
							}
							if (maskVal > 0)
								rhs(temp * d_in, ax, local_ind, output, param.noSensImage, SensImage, TOFWeights, param.TOF, param.nBins);
						}
						if (ii == iiEnd)
							break;
					}
					if (fp == 1) {
						if (nRays == 1) {
							for (size_t to = 0; to < param.nBins; to++) {
								forwardProjectAF(output, ax, ind, temp, to, CT);
								if (param.TOF)
									ind += nMeas;
							}
						}
						else {
							for (size_t to = 0; to < param.nBins; to++)
								ax[to + (size_t)param.nBins * (size_t)lor] *= temp;
						}
					}
				}
				else {
					int32_t tempi = 0, tempj = 0, tempk = 0;
					T txu = (T)0., tyu = (T)0., tzu = (T)0., tc = (T)0., tx0 = (T)1e8, ty0 = (T)1e8, tz0 = (T)1e8;
					// The ray parameter where the ray exits the FOV
					T tMax = (T)0.;
					bool skip = false;

					// Determine the above values and whether the ray intersects the FOV
					// Both detectors are on the same ring, but not perpendicular
					if (std::fabs(z_diff) < (T)1e-6) {
						tempk = tieFloor((detectors.zs - param.bz) / param.dz);
						if (tempk < 0 || tempk >= static_cast<int32_t>(param.Nz))
							continue;
						skip = siddon_pre_loop_2D_clamp(param.bx, param.by, x_diff, y_diff, bmaxx, bmaxy, param.dx, param.dy, param.Nx, param.Ny, tempi, tempj, txu, tyu, Np,
							detectors.ys, detectors.xs, detectors.yd, detectors.xd, tc, ux, uy, tx0, ty0, tMax);
					}
					//Detectors on different rings (e.g. oblique sinograms)
					else if (std::fabs(y_diff) < parTol) {
						tempj = perpendicular_start_tie(param.by, detectors.yd, param.dy, param.Ny);
						skip = siddon_pre_loop_2D_clamp(param.bx, param.bz, x_diff, z_diff, bmaxx, bmaxz, param.dx, param.dz, param.Nx, param.Nz, tempi, tempk, txu, tzu, Np,
							detectors.zs, detectors.xs, detectors.zd, detectors.xd, tc, ux, uz, tx0, tz0, tMax);
						if (detectors.yd > bmaxy || detectors.yd < param.by)
							skip = true;
					}
					else if (std::fabs(x_diff) < parTol) {
						tempi = perpendicular_start_tie(param.bx, detectors.xd, param.dx, param.Nx);
						skip = siddon_pre_loop_2D_clamp(param.by, param.bz, y_diff, z_diff, bmaxy, bmaxz, param.dy, param.dz, param.Ny, param.Nz, tempj, tempk, tyu, tzu, Np,
							detectors.zs, detectors.ys, detectors.zd, detectors.yd, tc, uy, uz, ty0, tz0, tMax);
						if (detectors.xd > bmaxx || detectors.xd < param.bx)
							skip = true;
					}
					else {
						skip = siddon_pre_loop_3D_clamp(param.bx, param.by, param.bz, x_diff, y_diff, z_diff, bmaxx, bmaxy, bmaxz, param.dx, param.dy, param.dz, param.Nx, param.Ny, param.Nz, tempi, tempj, tempk, txu, tyu, tzu,
							Np, detectors, tc, ux, uy, uz, tx0, ty0, tz0, tMax);
					}

					// Skip if the LOR does not intersect with the FOV
					if (skip) {
						continue;
					}
					// The ray is traversed from the source, the TOF coordinate starts at the source if the source is inside the FOV
					if (param.TOF) {
						TOFDis(std::max(tc, (T)0.), L, D);
						TOFInit(D, param.TOFCenters, invSigma, TOFF, param.nBins);
					}
					// The length of the ray inside the FOV (useTotLength == false, GPU: LL without TOTLENGTH)
					T LL = L;
					if (!param.useTotLength && !CT) {
						LL = L * (std::min(tMax, (T)1.) - std::max(tc, (T)0.));
						if (LL <= (T)0.)
							LL = L;
					}
						uint32_t local_ind = 0u;
						// The attenuation of the whole ray in the backprojection (GPU: ATN and BP). The same clipped traversal as below, i.e. only the part of
						// the ray between the end points t = [0, 1] is included
						if (param.attenuationCorrection && fp == 2 && param.CTAttenuation && !SPECT) {
							T tx0_c = tx0, ty0_c = ty0, tz0_c = tz0, tc_c = tc;
							int tempi_c = tempi, tempj_c = tempj, tempk_c = tempk;
							for (uint32_t ii = 0u; ii < Np; ii++) {
								local_ind = compute_ind(tempj_c, tempi_c, tempk_c, param.Nx, Nyx);
								bool passAtn = false;
								if (tz0_c < ty0_c && tz0_c < tx0_c)
									local_ele = compute_element_clip(tz0_c, tc_c, L, tzu, uz, tempk_c, passAtn);
								else if (ty0_c < tx0_c)
									local_ele = compute_element_clip(ty0_c, tc_c, L, tyu, uy, tempj_c, passAtn);
								else
									local_ele = compute_element_clip(tx0_c, tc_c, L, txu, ux, tempi_c, passAtn);
								compute_attenuation(local_ele, local_ind, param.atten, jelppi);
								// Unsigned comparisons, a negative index is a large unsigned number
								if (static_cast<uint32_t>(tempi_c) >= param.Nx || static_cast<uint32_t>(tempj_c) >= param.Ny || static_cast<uint32_t>(tempk_c) >= param.Nz)
									break;
							}
						}
					if (!CT) {
						if (param.projType == 1) {
							// Total length (useTotLength, for SPECT the length inside the support) or the length inside the FOV (GPU)
							temp = (T)1. / (LL * (T)(nRays));
						}
						if (param.attenuationCorrection && fp == 2 && param.CTAttenuation && !SPECT)
							temp *= std::exp(jelppi);
						// Measurement-based attenuation (GPU -DATNM), atten is then the attenuation of the measurements of the current subset
						else if (param.attenuationCorrection && !param.CTAttenuation)
							temp *= param.atten[lo];
						if (param.normalizationCorrection)
							temp *= local_norm;
						if (param.scatterCorrectionMult)
							temp *= local_scat;
						temp *= param.globalFactor;
					} else if (param.projType == 1 && nRays > 1) {
                        temp = (T)1. / (T)(nRays);
                    }
                    for (uint32_t ii = 0u; ii < Np; ii++) {
                        local_ind = compute_ind(tempj, tempi * d_N2, tempk, d_N3, Nyx);
                        localIndX = tempi;
                        localIndY = tempj;
                        localIndZ = tempk;
                        bool pass = false;
                        // The length of the intersection, clipped to the segment t = [0, 1] between the end points, and the move to the next voxel
                        // (GPU: compute_element)
                        if (tz0 < ty0 && tz0 < tx0)
                            local_ele = compute_element_clip(tz0, tc, L, tzu, uz, tempk, pass);
                        else if (ty0 < tx0)
                            local_ele = compute_element_clip(ty0, tc, L, tyu, uy, tempj, pass);
                        else
                            local_ele = compute_element_clip(tx0, tc, L, txu, ux, tempi, pass);
                        T local_ele2 = local_ele;
                        uint32_t local_ind2 = local_ind;
                        if (param.attenuationCorrection && (fp == 1 || SPECT) && param.CTAttenuation && pass) {
                            compute_attenuation(local_ele2, local_ind2, param.atten, jelppi);
							if (SPECT && param.projType == 1)
								local_ele *= std::exp(jelppi);
                        }
                        // Also moves D to the end of the voxel
                        if (param.TOF)
                            TOFLoop(local_ele2, param.TOFCenters, invSigma, D, TOFWeights, TOFF, param.nBins);
                        if (local_ele > (T)0.) {
                            if (fp == 1) {
                                denominator(ax, local_ind, local_ele, input, param.TOF, TOFWeights, param.nBins, lor, nRays);
                            }
                            else if (fp == 2) {
                                if (param.useMaskBP) {
                                    maskVal = param.maskBP[localIndX * d_N2 + localIndY * d_N3];
                                }
                                if (maskVal > 0)
                                    rhs(local_ele * temp, ax, local_ind, output, param.noSensImage, SensImage, TOFWeights, param.TOF, param.nBins);
                            }
                        }
                        // Unsigned comparisons, a negative index is a large unsigned number
                        if (static_cast<uint32_t>(tempi) >= param.Nx || static_cast<uint32_t>(tempj) >= param.Ny || static_cast<uint32_t>(tempk) >= param.Nz)
                            break;
                    }
                    if (param.attenuationCorrection && fp == 1 && param.CTAttenuation && !SPECT) {
                        temp *= std::exp(jelppi);
                    }
                    if (fp == 1) {
                        if (nRays == 1) {
                            for (size_t to = 0; to < param.nBins; to++) {
                                forwardProjectAF(output, ax, ind, temp, to, CT);
                                if (param.TOF)
                                    ind += nMeas;
                            }
                        }
                        else {
                            for (int64_t to = 0; to < param.nBins; to++)
                                ax[to + (int64_t)param.nBins * (int64_t)lor] *= temp;
                        }
                    }
                }
            }
        }
		if (nRays > 1 && fp == 1) {
			for (size_t to = 0; to < param.nBins; to++) {
				T apu = (T)0.;
				for (size_t kk = 0; kk < nRays; kk++) {
					apu += ax[to + param.nBins * kk];
				}
				ax[to] = apu;
			}
			for (size_t to = 0; to < param.nBins; to++) {
				forwardProjectAF(output, ax, ind, (T)1., to, CT);
				if (param.TOF)
					ind += nMeas;
			}
		}
	}
}
#ifdef _OPENMP
}
#endif

#ifndef AF
void improved_siddon_precomputed(paramStruct<double>& param, const int64_t nMeas, const double* x, const double* z, double* elements, size_t* indices, const uint16_t* lor1, const uint64_t* lor2,
	const bool CT = false, const uint16_t* detIndex = nullptr, const uint32_t nCores = 0);

void improved_siddon_precomputation_phase(paramStruct<double>& param, const int64_t nMeas, const double* x, const double* z, uint16_t* lor, const bool CT = false, 
	const uint16_t* detIndex = nullptr, const uint32_t nCores = 0);
#endif