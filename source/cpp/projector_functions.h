/**************************************************************************
* Header file for the various functions required by the improved Siddon, 
* orthogonal distance-based and volume-based ray tracers.
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

// Normalized distances below this are discarded in orthogonal ray tracer
#define H_THR 0.99
#define TOF_THR 0.0001
#define _2PI 0.3989422804014327
#ifndef TRAPZ_BINS
#define TRAPZ_BINS 4.
#endif
#define THR 0.01
#define CC 1e3

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
	// Projector type, 1 means (improved) Sidon, 2 is orthogonal and 3 the volume of intersection
	int projType = 1;
	// Radius of the tube of response with volume or width of the FWHM of the orthogonal (REQUIRED FOR ORTHOGONAL OR VOLUME)
	T orthWidth = static_cast<T>(1.);
	// Coordinates of the center of the voxels in the x-direction (REQUIRED FOR ORTHOGONAL OR VOLUME)
	T* x_center = nullptr;
	// Coordinates of the center of the voxels in the y-direction (REQUIRED FOR ORTHOGONAL OR VOLUME)
	T* y_center = nullptr;
	// Coordinates of the center of the voxels in the z-direction (REQUIRED FOR ORTHOGONAL OR VOLUME)
	T* z_center = nullptr;
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
	// The smallest orthogonal distance from the center of the ray to the center of the current voxel where the volume of intersection needs to be computed (REQUIRED FOR VOLUME)
	// if the volume is smaller than this value, the whole volume of the sphere is used
	T bmin = static_cast<T>(1.);
	// The largest allow orthogonal distance value. If the orthogonal distance from the ray to the center of the voxel is greater than this, the volume is not computed (REQUIRED FOR VOLUME)
	T bmax = static_cast<T>(1.);
	// The volume of the sphere (REQUIRED FOR VOLUME)
	T Vmax = static_cast<T>(1.);;
	// The precomputed volume of intersection values based on the orthogonal distance (REQUIRED FOR VOLUME)
	T* V = nullptr;
	// Are multiple layers used in PET
	bool nLayers = false;
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
	// Amount of forward projection masks
	int64_t numMaskFP = 1;
	// Total amount of projections (sinograms)
	int64_t nProjectionsGlobal = 1;
	// Collimator hole length (mm)
	//T colL = (T)1.;
	// Collimator hole diameter (detector end) [mm]
	//T collimatorBottomRadius = (T)1.;
	// Collimator hole diameter (source end) [mm]
	//T collimatorTopRadius = (T)1.;
	// Septal width (mm)
	//T dSeptal = (T)1.;
	// Number of rays traced per collimator hole
	//T nRaySPECT = (T)1.;
	// Collimator hexagon orientation: 1=vertical diameter smaller, 2=horizontal diameter smaller
	//T hexOrientation = (T)1.;
	// Method for tracing rays inside collimator hole: 1 for accurate location of rays, 2 for one cone at center of pixel
	//T coneMethod = (T)3.;
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
	if (subsetType >= 8 || subsets == 1 || subsetType == 1 || subsetType == 2 || subsetType == 4 || subsetType == 5) {
		const int64_t id = (ix + iy * size_x) * 4LL;
		if (nLayers) {
			const int64_t idz = iz * 3LL;
			const int layer = static_cast<int>(z[idz]);
			detectors.xs = x[id + layer * size_x * size_y];
			detectors.ys = x[id + 1 + layer * size_x * size_y];
			detectors.xd = x[id + 2 + layer * size_x * size_y];
			detectors.yd = x[id + 3 + layer * size_x * size_y];
			detectors.zs = z[idz + 1];
			detectors.zd = z[idz + 2];
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
			const int layer = static_cast<int>(z[indz]);
			detectors.xs = x[ind + layer * size_x * size_y];
			detectors.ys = x[ind + 1 + layer * size_x * size_y];
			detectors.xd = x[ind + 2 + layer * size_x * size_y];
			detectors.yd = x[ind + 3 + layer * size_x * size_y];
			detectors.zs = z[indz + 1];
			detectors.zd = z[indz + 2];
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

template <typename T>
inline void multirayCoordinateShiftXY(Det<T>& detectors, const int lor, const T cr, const uint16_t n_rays) {
	T interval = cr / (static_cast<T>(n_rays * 2));
	detectors.xs += (interval - cr / (T)2.);
	detectors.xd += (interval - cr / (T)2.);
	detectors.ys += (interval - cr / (T)2.);
	detectors.yd += (interval - cr / (T)2.);
	interval *= (T)2.;
	detectors.xs += interval * lor;
	detectors.xd += interval * lor;
	detectors.ys += interval * lor;
	detectors.yd += interval * lor;
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

template <typename T>
inline void get_detector_coordinates_CT(const T* x,const T* z, const uint32_t size_x, Det<T>& detectors, const int64_t lo, const uint32_t subsets,
	const uint32_t size_y, const int64_t ix, const int64_t iy, const int64_t iz, const T dPitch, const int64_t nProjections,
	const uint8_t list_mode_format, const bool pitch) {
	if (list_mode_format > 0) {
		detectors.xs = x[lo * 6];
		detectors.ys = x[lo * 6 + 1];
		detectors.zs = x[lo * 6 + 2];
		detectors.xd = x[lo * 6 + 3];
		detectors.yd = x[lo * 6 + 4];
		detectors.zd = x[lo * 6 + 5];
	}
	else {
		const int id = iz * 6;
		detectors.xs = x[id];
		detectors.ys = x[id + 1];
		detectors.zs = x[id + 2];
		detectors.xd = x[id + 3];
		detectors.yd = x[id + 4];
		detectors.zd = x[id + 5];
		const T indeksi1 = static_cast<T>(ix) - static_cast<T>(size_x) / 2. + .5;
		const T indeksi2 = static_cast<T>(iy) - static_cast<T>(size_y) / 2. + .5;
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

template <typename T>
inline void get_detector_coordinates_SPECT(const T* x, const T* z, Det<T>& detectors, const int64_t lo, const uint32_t lorXY, const uint32_t size_x, const uint32_t size_y, const int64_t ix, const int64_t iy, const int64_t iz, const T dPitchXY, const uint8_t list_mode_format, const uint32_t nRays2D, const T* rayShiftsDetector, const T* rayShiftsSource) {
	const int idx = iz * 6;
	detectors.xs = x[idx + 0]; // Detector panel center inner normal vector
	detectors.ys = x[idx + 1];
	detectors.zs = x[idx + 2];
	detectors.xd = x[idx + 3];
	detectors.yd = x[idx + 4];
	detectors.zd = x[idx + 5];

	const int idz = iz * 2; // Move ray to detector pixel
	const T indeksi1 = static_cast<T>(ix) - static_cast<T>(size_x) / 2. + .5;
	const T indeksi2 = static_cast<T>(iy) - static_cast<T>(size_y) / 2. + .5;
	detectors.xs += z[idz] * indeksi1;
	detectors.ys += z[idz + 1] * indeksi1;
	detectors.zs += dPitchXY * indeksi2;
	detectors.xd += z[idz] * indeksi1;
	detectors.yd += z[idz + 1] * indeksi1;
	detectors.zd += dPitchXY * indeksi2;

	/*if (lorXY == 0 && lo == 0 && false) {
		mexPrintf("dPitchXY = %f\n", (T)dPitchXY);
		mexPrintf("ix = %f\n", (T)ix);
		mexPrintf("iy = %f\n", (T)iy);
		mexPrintf("iz = %f\n", (T)iz);
		mexPrintf("indeksi1 = %f\n", (T)indeksi1);
		mexPrintf("indeksi2 = %f\n", (T)indeksi2);
		mexPrintf("detectors.xs = %f\n", (T)detectors.xs);
		mexPrintf("detectors.ys = %f\n", (T)detectors.ys);
		mexPrintf("detectors.zs = %f\n", (T)detectors.zs);
		mexPrintf("detectors.xd = %f\n", (T)detectors.xd);
		mexPrintf("detectors.yd = %f\n", (T)detectors.yd);
		mexPrintf("detectors.zd = %f\n", (T)detectors.zd);
	}*/
	
	if (nRays2D > 1) { // Add ray shift
		const int idr = lorXY * 2;

		detectors.xd += z[idz+0] * rayShiftsDetector[idr] / 2.;
		detectors.yd += z[idz+1] * rayShiftsDetector[idr] / 2.;
		detectors.zd += dPitchXY * rayShiftsDetector[idr + 1] / 2.;
		detectors.xs += z[idz+0] * rayShiftsSource[idr] / 2.;
		detectors.ys += z[idz+1] * rayShiftsSource[idr] / 2.;
		detectors.zs += dPitchXY * rayShiftsSource[idr + 1] / 2.;
	}	
	// Extend ray
	detectors.xs += 100 * (detectors.xs - detectors.xd);
	detectors.ys += 100 * (detectors.ys - detectors.yd);
	detectors.zs += 100 * (detectors.zs - detectors.zd);
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

template <typename T>
inline T normPDF(const T x, const T mu, const T sigma) {

	const T a = (x - mu) / sigma;

	return (static_cast<T>(_2PI) / sigma * std::exp((-0.5 * a * a)));
}

template <typename T>
inline T TOFWeight(const T element, const T sigma_x, const T D, const T DD, const T TOFCenter, T dX) {
	T output = normPDF(D, TOFCenter, sigma_x);
	dX = std::copysign(dX, DD);
	for (int tr = 1; tr < (int)(TRAPZ_BINS)-1; tr++)
		output += (normPDF(D - dX * (float)(tr), TOFCenter, sigma_x) * (T)2.);
	output += normPDF(D - std::copysign(element, DD), TOFCenter, sigma_x);
	return output;
}

template <typename T>
inline T TOFLoop(const T DD, const T element, const T* TOFCenter, const T sigma_x, T& D, const T epps, const uint32_t nBins) {
	T TOFSum = (T)0.;
	const T dX = element / (T)(TRAPZ_BINS - 1.);
	for (uint32_t to = 0; to < nBins; to++) {
		const T apu = TOFWeight(element, sigma_x, D, DD, TOFCenter[to], dX) * dX;
		TOFSum += apu;
	}
	if (TOFSum < epps)
		TOFSum = epps;
	return TOFSum;
}

template <typename T>
inline void TOFDis(const T x_diff, const T y_diff, const T z_diff, const T tc, const T LL, T& D, T& DD) {
	const T xI = x_diff * tc;
	const T yI = y_diff * tc;
	const T zI = z_diff * tc;
	D = norm(xI, yI, zI) - LL / (T)2.;
	DD = D;
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

template <typename T>
inline void denominator(std::vector<T>& ax, const uint32_t local_ind, T local_ele, const T* input, const bool TOF, const T element, const T TOFSum,
	const T DD, const T* TOFCenter, const T sigma_x, T& D, const uint32_t nBins, const int lor, const uint16_t nRays, const int projType) {
	T apu = (T)0.;
	forwardProject(local_ele, apu, local_ind, input);
	if (TOF) {
		const T dX = element / (T)(TRAPZ_BINS - 1.);
		for (uint32_t to = 0; to < nBins; to++) {
			const T joku = (TOFWeight(element, sigma_x, D, DD, TOFCenter[to], dX) * dX);
			if (nRays > 1)
				ax[(size_t)to + (size_t)nBins * (size_t)lor] += apu * joku / TOFSum;
			else
				ax[to] += apu * joku / TOFSum;
		}
	}
	else {
		if (nRays > 1)
			ax[lor] += apu;
		else
			ax[0] += apu;
	}
}

// Compute the backprojection
template <typename T>
inline void rhs(const T local_ele, const std::vector<T>& ax, const uint32_t local_ind, T* output, const bool no_norm, T* sensImage, const T element,
	const T sigma_x, T& D, const T DD, const T* TOFCenter, const T TOFSum, const bool TOF, const uint32_t nBins, const int projType) {
	T yaxTOF = (T)0.;
	T val = (T)0.;
	if (TOF) {
		const T dX = element / (TRAPZ_BINS - (T)1.);
		for (uint32_t to = 0; to < nBins; to++) {

			const T apu = local_ele * ((TOFWeight(element, sigma_x, D, DD, TOFCenter[to], dX) * dX) / TOFSum);

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
	const int TYPE, const T ys, const T xs, const T yd, const T xd, T& tc, int32_t& u1, int32_t& u2, T& t10, T& t20, const int projType = 1, bool& xy = false) {
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
	if (projType > 1) {
		if (tc == t10 || tc == txback)
			xy = true;
		else
			xy = false;
	}

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
	uint32_t& Np, const int TYPE, const Det<T> detectors, T& tc, int32_t& iu, int32_t& ju, int32_t& ku, T& tx0, T& ty0, T& tz0, const int projType = 1, bool& xy = false) {

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
	if (projType > 1) {
		if (tc == tx0 || tc == txback)
			xy = true;
		else
			xy = false;
	}

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

// Compute the orthogonal distance from the ray to the current voxel (center)
// For orthogonal distance-based ray tracer, the distance is normalized
template <typename T>
inline T compute_element_orth_3D(const T xs, const T ys, const T zs, const T xl, const T yl, const T zl, const T crystal_size_z,
	const T xp, const int projType) {

	const T x0 = xp - xs;

	// Cross product
	const T y1 = zl * x0 - xl;
	const T z1 = -yl * x0 + ys;

	const T normi = norm(zs, y1, z1);

	if (projType == 3)
		return (normi / crystal_size_z);
	else
		return (1.f - normi / crystal_size_z);
}

// compute voxel index, orthogonal distance based or volume of intersection ray tracer
inline uint32_t compute_ind_orth_3D(const uint32_t tempi, const uint32_t tempijk, const int tempk, const uint32_t d_N, const uint32_t Nyx) {
	uint32_t local_ind = tempi * d_N + tempijk + (uint32_t)(tempk)*Nyx;
	return local_ind;
}

// This function computes either the forward or backward projection for the current voxel
// The normalized orthogonal distance or the volume of the (spherical) voxel is computed before the forward or backward projections
template <typename T>
inline bool orthogonalHelper3D(const uint32_t tempi, const int uu, const uint32_t d_N2, const uint32_t d_N3, const uint32_t d_Nxy, const int zz, const T s2, const T l3, const T l1, const T l2,
	const T diff1, const T diffZ, const T kerroin, const T center2, const T bmin, const T bmax, const T Vmax, T* V, const bool XY, std::vector<T>& ax, const T temp, const T* input,
	T* d_Summ, T* d_output, const bool no_norm, const T element, const T sigma_x, T& D, const T DD, const T* TOFCenter, const T TOFSum, const bool TOF, const uint8_t fp, const int projType,
	const uint32_t nBins, const int lor, const uint16_t nRays, const bool useMaskBP = false, const uint8_t* maskBP = nullptr, const T attApu = (T)0.f, const bool SPECT = false, 
	const bool attenuationCorrection = false) {
	T local_ele = compute_element_orth_3D(s2, l3, l1, l2, diff1, diffZ, kerroin, center2, projType);
	uint8_t maskVal = 1;
	if (projType == 3) {
		if (local_ele >= bmax) {
			return true;
		}
		if (local_ele < bmin)
			local_ele = Vmax;
		else
			local_ele = V[(uint32_t)(std::round((local_ele - bmin) * CC))];
	}
	else {
		if (local_ele <= THR) {
			return true;
		}
	}
	if (SPECT && attenuationCorrection)
		local_ele *= attApu;
	uint32_t local_ind = compute_ind_orth_3D(tempi, uu * d_N3, (zz), d_N2, d_Nxy);
	if (fp == 1) {
		denominator(ax, local_ind, local_ele, input, TOF, element, TOFSum, DD, TOFCenter, sigma_x, D, nBins, lor, nRays, projType);
	}
	else if (fp == 2) {
		if (useMaskBP)
			maskVal = maskBP[tempi + uu * d_N3];
		if (maskVal > 0)
			rhs(local_ele * temp, ax, local_ind, d_output, no_norm, d_Summ, element, sigma_x, D, DD, TOFCenter, TOFSum, TOF, nBins, projType);
	}
	return false;
}

// Both the orthogonal and volumme of intersection ray tracers loop through all the neighboring voxels of the current voxel
// Both also proceed through each X or Y slice, depending on the incident direction
// This function simply loops through each X or Y voxel and Z voxels in the current slice
// Forward or backward projection is computed in the helper function
template <typename T>
inline int orthDistance3D(const uint32_t tempi, const T diff1, const T diff2, const T diffZ, const T center1, const T* center2, const T* centerZ, const T temp, int temp2, const int tempk,
	const T s1, const T s2, const T sZ, const uint32_t d_Nxy, const T kerroin, const uint32_t d_N1, const uint32_t d_N2, const uint32_t d_N3, const uint32_t d_Nz, const T bmin,
	const T bmax, const T Vmax, T* V, const bool XY, std::vector<T>& ax, const T* input, const bool no_norm, T* Summ, T* output, const T element, const T sigma_x, T& D, const T DD,
	T* TOFCenter, const T TOFSum, const bool TOF, const uint8_t fp, const int projType, const uint32_t nBins, const int lor, const uint16_t nRays, const bool useMaskBP = false, const uint8_t* maskBP = nullptr, 
	const T attApu = (T)0.f, const bool SPECT = false, const bool attenuationCorrection = false) {
	int uu = 0;
	bool breikki = false;
	// y0
	const T v0 = center1 - s1;
	// xl * y0
	const T l3 = diff1 * v0;
	// zl * y0
	const T apu1 = diffZ * v0;
	const int maksimiZ = (int)(d_Nz);
	const int minimiZ = 0;
	const int maksimiXY = (int)(d_N1);
	const int minimiXY = 0;
	int uu1 = 0, uu2 = 0;
	for (int zz = tempk; zz < maksimiZ; zz++) {
		// z0
		const T z0 = centerZ[zz] - sZ;
		// x1 = yl * z0 - zl * y0
		const T l1 = diff2 * z0 - apu1;
		// xl * z0
		const T l2 = diff1 * z0;
		for (uu1 = temp2; uu1 < maksimiXY; uu1++) {
			breikki = orthogonalHelper3D(tempi, uu1, d_N2, d_N3, d_Nxy, zz, s2, l3, l1, l2, diff2, diffZ, kerroin, center2[uu1], bmin, bmax, Vmax, V,
				XY, ax, temp, input, Summ, output, no_norm, element, sigma_x, D, DD, TOFCenter, TOFSum, TOF, fp, projType, nBins, lor, nRays, useMaskBP, 
				maskBP, attApu, SPECT, attenuationCorrection);
			if (breikki) {
				break;
			}
			uu++;
		}
		for (uu2 = temp2 - 1; uu2 >= minimiXY; uu2--) {
			breikki = orthogonalHelper3D(tempi, uu2, d_N2, d_N3, d_Nxy, zz, s2, l3, l1, l2, diff2, diffZ, kerroin, center2[uu2], bmin, bmax, Vmax, V,
				XY, ax, temp, input, Summ, output, no_norm, element, sigma_x, D, DD, TOFCenter, TOFSum, TOF, fp, projType, nBins, lor, nRays, useMaskBP, 
				maskBP, attApu), SPECT, attenuationCorrection;
			if (breikki) {
				break;
			}
			uu++;
		}
		if (uu1 == temp2 && uu2 == temp2 - 1 && breikki)
			break;
	}
	for (int zz = tempk - 1; zz >= minimiZ; zz--) {
		const T z0 = centerZ[zz] - sZ;
		const T l1 = diff2 * z0 - apu1;
		const T l2 = diff1 * z0;
		for (uu1 = temp2; uu1 < maksimiXY; uu1++) {
			breikki = orthogonalHelper3D(tempi, uu1, d_N2, d_N3, d_Nxy, zz, s2, l3, l1, l2, diff2, diffZ, kerroin, center2[uu1], bmin, bmax, Vmax, V,
				XY, ax, temp, input, Summ, output, no_norm, element, sigma_x, D, DD, TOFCenter, TOFSum, TOF, fp, projType, nBins, lor, nRays, useMaskBP, 
				maskBP, attApu, SPECT, attenuationCorrection);
			if (breikki) {
				break;
			}
			uu++;
		}
		for (uu2 = temp2 - 1; uu2 >= minimiXY; uu2--) {
			breikki = orthogonalHelper3D(tempi, uu2, d_N2, d_N3, d_Nxy, zz, s2, l3, l1, l2, diff2, diffZ, kerroin, center2[uu2], bmin, bmax, Vmax, V,
				XY, ax, temp, input, Summ, output, no_norm, element, sigma_x, D, DD, TOFCenter, TOFSum, TOF, fp, projType, nBins, lor, nRays, useMaskBP, 
				maskBP, attApu, SPECT, attenuationCorrection);
			if (breikki) {
				break;
			}
			uu++;
		}
		if (uu1 == temp2 && uu2 == temp2 - 1 && breikki)
			break;
	}
	return uu;
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

	// if 0, then determines whether the LOR intercepts the FOV (i.e. no precomputation phase performed)
	const static int TYPE = 0;

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
//#ifdef _OPENMP
//	size_t threads = omp_get_max_threads();
//	if (DEBUG)
//		mexPrintf("threadsOMP = %u\n", threads);
//#else
//	size_t threads = 1ULL;
//	if (DEBUG)
//		mexPrintf("threads = %u\n", threads);
//#endif
	//if (DEBUG)
	//	mexPrintf("nMeas = %u\n", nMeas);

	uint32_t nRays = param.nRays2D * param.nRays3D;
#ifdef _OPENMP
#pragma omp parallel
	{
		std::vector<T> ax(param.nBins);
#if _OPENMP >= 201511 && defined(MATLAB)
#pragma omp for schedule(monotonic:dynamic, nChunks)
#else
#pragma omp for schedule(dynamic, nChunks)
#endif
#else
	std::vector<T> ax(param.nBins);
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
			/*if (lo == 0) {
				mexPrintf("ix = %f\n", (T)ix);
				mexPrintf("iy = %f\n", (T)iy);
				mexPrintf("iz = %f\n", (T)iz);
			}*/
		}
		else {
			if ((param.listMode > 0 && param.computeSensIm == 0) || (!(CT || SPECT) && param.listMode == 0)) {
				ix = lo;
				iy = 0;
				iz = 0;
			}
			else if (param.listMode > 0 && param.computeSensIm) {
				ix = lo % param.det_per_ring;
				iy = (lo / param.det_per_ring) % param.det_per_ring;
				iz = lo / (param.det_per_ring * param.det_per_ring) % param.rings;
				izD = lo / (param.det_per_ring * param.det_per_ring) / param.rings;
			}
			else {
				iz = ix / ((int64_t)param.size_x * (int64_t)param.size_y);
				iy = ix % param.size_y;
				ix = (ix / param.size_y) % param.size_x;
			}
		}
		int64_t ind = lo;
		if (param.useMaskFP) {
			const int64_t idx1 = ix + iy * (int64_t)param.size_x;
			int64_t idx2 = 0;
			if (param.numMaskFP > 1) { // Shift to the correct detector panel mask 
				int64_t numProjPerDetector = param.nProjectionsGlobal / param.numMaskFP;
				int64_t currentDetector = iz / numProjPerDetector;
				idx2 = (int64_t)param.size_x * (int64_t)param.size_x * currentDetector;
			}
			const bool maskVal = param.maskFP[idx1+idx2];
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
		T local_norm = (T)0.;
		T local_scat = (T)0.;
		if (param.normalizationCorrection)
			local_norm = static_cast<T>(param.normCoefs[lo]);
		if (param.scatterCorrectionMult)
			local_scat = static_cast<T>(param.scatterCoefs[lo]);
		int lor = -1;

		// Loop through the rays
		//for (int currentShift = 0; currentShift < nShift; currentShift++) {
		for (int lorZ = 0; lorZ < param.nRays3D; lorZ++) {
			for (int lorXY = 0; lorXY < param.nRays2D; lorXY++) {
				lor++;

				//if (!SPECT) {
				//Det<T> detectors;
				detectors.xd = (T)0., detectors.xs = (T)0., detectors.yd = (T)0., detectors.ys = (T)0., detectors.zd = (T)0., detectors.zs = (T)0.;
				if (CT) { // CT data
					get_detector_coordinates_CT(x, z, param.size_x, detectors, lo, param.subsets, param.size_y, ix, iy, iz, param.dPitchZ, param.nProjections, param.listMode, param.pitch);
				} else if (SPECT) { // SPECT data
					get_detector_coordinates_SPECT(x, z, detectors, lo, lorXY, param.size_x, param.size_y, ix, iy, iz, param.dPitchXY, param.listMode, param.nRays2D, param.rayShiftsDetector, param.rayShiftsSource);	
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
				if ((y_diff == 0. && x_diff == 0. && z_diff == 0.) || (y_diff == 0. && x_diff == 0.) || std::isinf(y_diff) || std::isinf(x_diff)) {
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
					uint32_t local_ind = 0u;
					int localIndX = 0, localIndY = 0, localIndZ = 0;
					T D = 0., DD = 0.;
					T local_ele = 0.;
					bool XY = false;
					T kerroin = 0.;
					T TotV = 0.;
					T dI = 0., TOFSum = 0.;
					T* center2 = nullptr, * center1 = nullptr;
					uint8_t maskVal = 1;
					if (param.projType == 2)
						kerroin = L * param.orthWidth;
					else if (param.projType == 3) {
						kerroin = L;
						TotV = L * (T)(M_PI) * param.orthWidth * param.orthWidth;
					}

				if (std::fabs(z_diff) < 1e-8 && (std::fabs(y_diff) < 1e-8 || std::fabs(x_diff) < 1e-8)) {

					// Ring number
					int32_t tempk = static_cast<int32_t>(std::fabs(detectors.zs - param.bz) / param.dz);
					if (tempk < 0 || tempk >= param.Nz)
						continue;
					int indO = 0;
					T d_b, dd, d_db, d_d2;

					// Detectors are perpendicular
					// Siddon cannot be applied --> trivial to compute
					int32_t apuX1, apuX2;
					T dT1, dT2;
					if (std::fabs(y_diff) < 1e-8 && detectors.yd <= bmaxy && detectors.yd >= param.by) {
						apuX1 = 0;
						apuX2 = param.Nx - 1;
						dT1 = param.dx;
						dT2 = param.dx;
						d_b = param.by;
						dd = detectors.yd;
						d_db = param.dy;
						if (param.projType > 1) {
							center2 = param.x_center;
							center1 = param.y_center;
						}
						if (param.projType == 1) {
							T dist1, dist2 = 0.f;
							if (detectors.xs > detectors.xd) {
								dist1 = (param.bx - detectors.xd);
								dist2 = (param.bx + static_cast<T>(param.Nx) * param.dx - detectors.xs);
							}
							else {
								dist1 = (param.bx - detectors.xs);
								dist2 = (param.bx + static_cast<T>(param.Nx) * param.dx - detectors.xd);
							}
							for (int kk = 0; kk < param.Nx; kk++) {
								if (dist1 >= (T)0.) {
									apuX1 = kk;
									if (kk == 0)
										dT1 = param.dx;
									else
										dT1 = std::min(dist1, param.dx);
									break;
								}
								dist1 += param.dx;
							}
							for (int kk = param.Nx - 1; kk >= apuX1; kk--) {
								if (dist2 <= (T)0.) {
									apuX2 = kk;
									if (kk == param.Nx - 1)
										dT2 = param.dx;
									else
										dT2 = std::min(-dist2, param.dx);
									break;
								}
								dist2 -= param.dx;
							}
						}
						XY = true;
						d_d2 = param.dx;
						d_N0 = param.Ny;
						d_N1 = param.Nx;
						d_N2 = param.Ny;
						d_N3 = 1u;
					}
					else if (std::fabs(x_diff) < 1e-8 && detectors.xd <= bmaxx && detectors.xd >= param.bx) {
						apuX1 = 0;
						apuX2 = param.Ny - 1;
						dT1 = param.dx;
						dT2 = param.dx;
						d_b = param.bx;
						dd = detectors.xd;
						if (param.projType > 1) {
							center2 = param.y_center;
							center1 = param.x_center;
						}
						if (param.projType == 1) {
							T dist1, dist2 = 0.f;
							if (detectors.ys > detectors.yd) {
								dist1 = (param.by - detectors.yd);
								dist2 = (param.by + static_cast<T>(param.Ny) * param.dy - detectors.ys);
							}
							else {
								dist1 = (param.by - detectors.ys);
								dist2 = (param.by + static_cast<T>(param.Ny) * param.dy - detectors.yd);
							}
							for (int kk = 0; kk < param.Ny; kk++) {
								if (dist1 >= (T)0.) {
									apuX1 = kk;
									if (kk == 0)
										dT1 = param.dy;
									else
										dT1 = std::min(dist1, param.dy);
									break;
								}
								dist1 += param.dy;
							}
							for (int kk = param.Ny - 1; kk >= apuX1; kk--) {
								if (dist2 <= (T)0.) {
									apuX2 = kk;
									if (kk == param.Ny - 1)
										dT2 = param.dy;
									else
										dT2 = std::min(-dist2, param.dy);
									break;
								}
								dist2 -= param.dy;
							}
						}
						d_d2 = param.dy;
						d_db = param.dx;
						T xs_apu = detectors.xs;
						detectors.xs = detectors.ys;
						detectors.ys = xs_apu;
						T xdiff_apu = x_diff;
						x_diff = y_diff;
						y_diff = xdiff_apu;
					}
					else
						continue;
					localIndZ = tempk;
					temp = perpendicular_elements(d_N2, dd, d_d2, d_b, d_db, d_N0, d_N1, param.atten, local_norm, param.attenuationCorrection, param.normalizationCorrection, 
						param.CTAttenuation, tempk, d_N3, param.globalFactor, param.scatterCorrectionMult, local_scat, localIndX, localIndY, localIndZ, L, nRays, 
						param.projType, param.Nx, param.Ny, CT, lo, SPECT);
					local_ind = tempk;
					if (param.projType == 3)
						temp *= ((T)1. / TotV);
					if (param.projType > 1 || (fp == 2 && param.useMaskBP)) {
						if (d_N2 == 1)
							indO = localIndX;
						else
							indO = localIndY;
					}
					if (param.TOF) {
						dI = (d_d2 * d_N1) / std::copysign(2.f, y_diff);
						D = dI;
						DD = D;
					}
					T attApu = (T)0.;

					for (uint32_t ii = apuX1; ii < apuX2; ii++) {
						T d_in = d_d2;
						if (ii == apuX1) {
							local_ind += d_N3 * ii;
						}
						if (apuX1 > 0 && ii == apuX1)
							d_in = dT1;
						else if (apuX2 < d_N1 - 1 && ii == apuX2)
							d_in = dT2;
						if (param.attenuationCorrection && SPECT && param.CTAttenuation) {
							compute_attenuation(d_in, local_ind, param.atten, jelppi);
							if (param.projType == 1)
								d_in *= std::exp(jelppi);
							else
								attApu = std::exp(jelppi);
						}
						if (param.TOF)
							TOFSum = TOFLoop(DD, d_d2, param.TOFCenters, param.sigma_x, D, param.epps, param.nBins);
						if (param.projType > 1) {
							orthDistance3D(ii, y_diff, x_diff, z_diff, center1[ii], center2, param.z_center, temp, indO, localIndZ, detectors.xs, detectors.ys, detectors.zs, Nyx, kerroin, d_N1, d_N3, d_N2, param.Nz, 
								param.bmin, param.bmax, param.Vmax, param.V, XY, ax, input, param.noSensImage, SensImage, output, d_d2, param.sigma_x, D, DD, param.TOFCenters, TOFSum, param.TOF, fp, param.projType, 
								param.nBins, lor, nRays, param.useMaskBP, param.maskBP, attApu, SPECT, param.attenuationCorrection);
						}
						else {
							if (fp == 1) {
								denominator(ax, local_ind, d_in, input, param.TOF, d_in, TOFSum, DD, param.TOFCenters, param.sigma_x, D, param.nBins, lor, nRays, param.projType);
							}
							else if (fp == 2) {
								if (param.useMaskBP) {
									maskVal = param.maskBP[indO * d_N2 + ii * d_N3];
								}
								if (maskVal > 0)
									rhs(temp * d_in, ax, local_ind, output, param.noSensImage, SensImage, d_in, param.sigma_x, D, DD, param.TOFCenters, TOFSum, param.TOF, param.nBins, param.projType);
							}
						}
						local_ind += d_N3;
						if (param.TOF)
							D -= std::copysign(d_d2, DD);
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
					bool skip = false, XY = true;

					// Determine the above values and whether the ray intersects the FOV
					// Both detectors are on the same ring, but not perpendicular
					if (std::fabs(z_diff) < (T)1e-8) {
						tempk = static_cast<int>(fabs(detectors.zs - param.bz) / param.dz);
						if (tempk < 0 || tempk >= param.Nz)
							continue;
						skip = siddon_pre_loop_2D(param.bx, param.by, x_diff, y_diff, bmaxx, bmaxy, param.dx, param.dy, param.Nx, param.Ny, tempi, tempj, txu, tyu, Np, TYPE,
							detectors.ys, detectors.xs, detectors.yd, detectors.xd, tc, ux, uy, tx0, ty0, param.projType, XY);
					}
					//Detectors on different rings (e.g. oblique sinograms)
					else if (std::fabs(y_diff) < (T)1e-8) {
						tempj = perpendicular_start(param.by, detectors.yd, param.dy, param.Ny);
						skip = siddon_pre_loop_2D(param.bx, param.bz, x_diff, z_diff, bmaxx, bmaxz, param.dx, param.dz, param.Nx, param.Nz, tempi, tempk, txu, tzu, Np, TYPE,
							detectors.zs, detectors.xs, detectors.zd, detectors.xd, tc, ux, uz, tx0, tz0, param.projType, XY);
						XY = true;
						if (detectors.yd > bmaxy || detectors.yd < param.by)
							skip = true;
					}
					else if (std::fabs(x_diff) < (T)1e-8) {
						tempi = perpendicular_start(param.bx, detectors.xd, param.dx, param.Nx);
						skip = siddon_pre_loop_2D(param.by, param.bz, y_diff, z_diff, bmaxy, bmaxz, param.dy, param.dz, param.Ny, param.Nz, tempj, tempk, tyu, tzu, Np, TYPE,
							detectors.zs, detectors.ys, detectors.zd, detectors.yd, tc, uy, uz, ty0, tz0, param.projType, XY);
						XY = false;
						if (detectors.xd > bmaxx || detectors.xd < param.bx)
							skip = true;
					}
					else {
						skip = siddon_pre_loop_3D(param.bx, param.by, param.bz, x_diff, y_diff, z_diff, bmaxx, bmaxy, bmaxz, param.dx, param.dy, param.dz, param.Nx, param.Ny, param.Nz, tempi, tempj, tempk, tyu, txu, tzu,
							Np, TYPE, detectors, tc, ux, uy, uz, tx0, ty0, tz0, param.projType, XY);
					}

					// Skip if the LOR does not intersect with the FOV
					if (skip) {
						continue;
					}
					if (param.TOF)
						TOFDis(x_diff, y_diff, z_diff, tc, L, D, DD);
					int tempi_b = 0, u_b = 0;
					uint32_t d_Nb = 0U;
					uint32_t local_ind = 0u;
					T t0_b = (T)0., tu_b = (T)0., diff_b = (T)0., xs = (T)0., ys = (T)0.;
					T attApu = (T)0.;
					uint32_t d_NNx = param.Nx;
					uint32_t d_NNy = param.Ny;
					uint32_t d_NNz = param.Nz;
					T tx0_c = tx0, ty0_c = ty0, tz0_c = tz0, txu_c = txu, tyu_c = tyu, tzu_c = tzu, tc_c = tc;
					int tempi_c = tempi, tempj_c = tempj, tempk_c = tempk, ux_c = ux, uy_c = uy, uz_c = uz;

					if (param.attenuationCorrection && fp == 2 && param.CTAttenuation && !SPECT) {
						T tc_a = tc;
						for (uint32_t ii = 0u; ii < Np; ii++) {
							local_ind = compute_ind(tempj_c, tempi_c, tempk_c, param.Nx, Nyx);
							if (tz0_c < ty0_c && tz0_c < tx0_c) {
								local_ele = compute_element(tz0_c, tc, L, tzu_c, uz_c, tempk_c);
							}
							else if (ty0_c < tx0_c) {
								local_ele = compute_element(ty0_c, tc, L, tyu_c, uy_c, tempj_c);
							}
							else {
								local_ele = compute_element(tx0_c, tc, L, txu_c, ux_c, tempi_c);
							}
							compute_attenuation(local_ele, local_ind, param.atten, jelppi);
							if (tempi_c < 0 || tempi_c >= param.Nx || tempj_c < 0 || tempj_c >= param.Ny || tempk_c < 0 || tempk_c >= d_NNz) {
								break;
							}
						}
						tc = tc_a;
						tx0_c = tx0, ty0_c = ty0, tz0_c = tz0, txu_c = txu, tyu_c = tyu, tzu_c = tzu, tc_c = tc;
						tempi_c = tempi, tempj_c = tempj, tempk_c = tempk, ux_c = ux, uy_c = uy, uz_c = uz;
					}
					if (param.projType > 1) {
						if (!XY) {
							tempi_b = tempi;
							tempi = tempj;
							tempj = tempi_b;
							xs = detectors.ys;
							ys = detectors.xs;
							u_b = ux;
							ux = uy;
							uy = u_b;
							t0_b = tx0;
							tx0 = ty0;
							ty0 = t0_b;
							tu_b = txu;
							txu = tyu;
							tyu = tu_b;
							diff_b = x_diff;
							x_diff = y_diff;
							y_diff = diff_b;
							center1 = param.y_center;
							center2 = param.x_center;
							d_NNx = param.Ny;
							d_NNy = param.Nx;
							d_N0 = param.Ny;
							d_N1 = param.Nx;
							d_N2 = param.Nx;
							d_N3 = 1;
						}
						else {
							xs = detectors.xs;
							ys = detectors.ys;
							center1 = param.x_center;
							center2 = param.y_center;
							d_N0 = param.Nx;
							d_N1 = param.Ny;
							d_N2 = 1;
							d_N3 = param.Nx;
						}
					}
					T tx0_a = tx0, ty0_a = ty0, tz0_a = tz0;
					int tempi_a = tempi, tempj_a = tempj, tempk_a = tempk;
					if (!CT) {
						if (param.projType == 3) {
							temp = (T)1. / TotV;
                        } else if (param.projType == 1) {
							temp = (T)1. / (L * (T)(nRays));
                        }
//						else if (param.projType == 1)
//							temp = (T)1. / L;
						if (param.attenuationCorrection && fp == 2 && param.CTAttenuation && !SPECT)
							temp *= std::exp(jelppi);
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
                    T L_SPECT = 0.f;
                    if (SPECT && (param.projType == 1)) {
                        temp *= (L * (T)(nRays));
                    }

                    for (uint32_t ii = 0u; ii < Np; ii++) {
                        local_ele = (T)0.;
                        local_ind = compute_ind(tempj, tempi * d_N2, tempk, d_N3, Nyx);
                        localIndX = tempi;
                        localIndY = tempj;
                        localIndZ = tempk;
                        if (param.projType > 1) {
                            tx0_a = tx0;
                            ty0_a = ty0;
                            tz0_a = tz0;
                            tempi_a = tempi;
                            tempj_a = tempj;
                            tempk_a = tempk;
                        }
						bool pass = false;
                        if (tz0 < ty0 && tz0 < tx0) {
                            if (tz0 >= (T)0. && tz0 <= (T)1.) {
                                if (tc < (T)0.) {
                                    local_ele = tz0 * L;
                                    compute_element(tz0, tc, L, tzu, uz, tempk);
                                }
                                else
                                    local_ele = compute_element(tz0, tc, L, tzu, uz, tempk);
								pass = true;
                            }
                            else if (tc >= (T)0. && tc <= (T)1.) {
                                if (tz0 > (T)1.) {
                                    local_ele = ((T)1. - tc) * L;
                                    compute_element(tz0, tc, L, tzu, uz, tempk);
                                }
                                else
                                    local_ele = compute_element(tz0, tc, L, tzu, uz, tempk);
								pass = true;
                            }
                            else
                                compute_element(tz0, tc, L, tzu, uz, tempk);
                        }
                        else if (ty0 < tx0) {
                            if (ty0 >= (T)0. && ty0 <= (T)1.) {
                                if (tc < (T)0.) {
                                    local_ele = ty0 * L;
                                    compute_element(ty0, tc, L, tyu, uy, tempj);
                                }
                                else
                                    local_ele = compute_element(ty0, tc, L, tyu, uy, tempj);
								pass = true;
                            }
                            else if (tc >= (T)0. && tc <= (T)1.) {
                                if (ty0 > (T)1.) {
                                    local_ele = ((T)1. - tc) * L;
                                    compute_element(ty0, tc, L, tyu, uy, tempj);
                                }
                                else
                                    local_ele = compute_element(ty0, tc, L, tyu, uy, tempj);
								pass = true;
                            }
                            else
                                compute_element(ty0, tc, L, tyu, uy, tempj);
                        }
                        else {
                            if (tx0 >= (T)0. && tx0 <= (T)1.) {
                                if (tc < (T)0.) {
                                    local_ele = tx0 * L;
                                    compute_element(tx0, tc, L, txu, ux, tempi);
                                }
                                else
                                    local_ele = compute_element(tx0, tc, L, txu, ux, tempi);
								pass = true;
                            }
                            else if (tc >= (T)0. && tc <= (T)1.) {
                                if (tx0 > (T)1.) {
                                    local_ele = ((T)1. - tc) * L;
                                    compute_element(tx0, tc, L, txu, ux, tempi);
                                }
                                else
                                    local_ele = compute_element(tx0, tc, L, txu, ux, tempi);
								pass = true;
                            }
                            else
                                compute_element(tx0, tc, L, txu, ux, tempi);
                        }
                        T local_ele2 = local_ele;
                        uint32_t local_ind2 = local_ind;
                        if (param.projType > 1 && ((param.attenuationCorrection && fp == 1 && param.CTAttenuation) || param.TOF)) {
                            if (param.attenuationCorrection && fp == 1 && param.CTAttenuation)
                                local_ind2 = compute_ind(tempj_c, tempi_c, tempk_c, param.Nx, Nyx);
                            if (tz0_c < ty0_c && tz0_c < tx0_c) {
                                local_ele2 = compute_element(tz0_c, tc_c, L, tzu_c, uz_c, tempk_c);
                            }
                            else if (ty0_c < tx0_c) {
                                local_ele2 = compute_element(ty0_c, tc_c, L, tyu_c, uy_c, tempj_c);
                            }
                            else {
                                local_ele2 = compute_element(tx0_c, tc_c, L, txu_c, ux_c, tempi_c);
                            }
                        }
                        if (param.attenuationCorrection && (fp == 1 || SPECT) && param.CTAttenuation && pass) {
                            compute_attenuation(local_ele2, local_ind2, param.atten, jelppi);
							if (SPECT && param.projType == 1)
								local_ele *= std::exp(jelppi);
							else if (SPECT)
								attApu = std::exp(jelppi);
                        }
                        if (param.TOF)
                            TOFSum = TOFLoop(DD, local_ele2, param.TOFCenters, param.sigma_x, D, param.epps, param.nBins);
                        if (param.projType > 1) {
                            if (ii == 0) {
                                if (ux >= 0) {
                                    for (int kk = tempi_a - 1; kk >= 0; kk--) {
                                        int uu = orthDistance3D(kk, y_diff, x_diff, z_diff, center1[kk], center2, param.z_center, temp, tempj_a, tempk_a, xs, ys, detectors.zs, Nyx, kerroin, d_N1, d_N2, d_N3,
                                            param.Nz, param.bmin, param.bmax, param.Vmax, param.V, XY, ax, input, param.noSensImage, SensImage, output, local_ele2, param.sigma_x, D, DD, param.TOFCenters, TOFSum, param.TOF, fp,
                                            param.projType, param.nBins, lor, nRays, param.useMaskBP, param.maskBP, attApu, SPECT, param.attenuationCorrection);
                                        if (uu == 0)
                                            break;
                                    }
                                }
                                else {
                                    for (int kk = tempi_a + 1; kk < d_NNx; kk++) {
                                        int uu = orthDistance3D(kk, y_diff, x_diff, z_diff, center1[kk], center2, param.z_center, temp, tempj_a, tempk_a, xs, ys, detectors.zs, Nyx, kerroin, d_N1, d_N2, d_N3,
                                            param.Nz, param.bmin, param.bmax, param.Vmax, param.V, XY, ax, input, param.noSensImage, SensImage, output, local_ele2, param.sigma_x, D, DD, param.TOFCenters, TOFSum, param.TOF, fp,
                                            param.projType, param.nBins, lor, nRays, param.useMaskBP, param.maskBP, attApu, SPECT, param.attenuationCorrection);
                                        if (uu == 0)
                                            break;
                                    }
                                }
                            }
                            if (tz0_a >= tx0_a && ty0_a >= tx0_a) {
                                orthDistance3D(localIndX, y_diff, x_diff, z_diff, center1[localIndX], center2, param.z_center, temp, localIndY, localIndZ, xs, ys, detectors.zs, Nyx, kerroin, d_N1, d_N2, d_N3,
                                    param.Nz, param.bmin, param.bmax, param.Vmax, param.V, XY, ax, input, param.noSensImage, SensImage, output, local_ele2, param.sigma_x, D, DD, param.TOFCenters, TOFSum, param.TOF, fp,
                                    param.projType, param.nBins, lor, nRays, param.useMaskBP, param.maskBP, attApu, SPECT, param.attenuationCorrection);
                            }
                        }
                        else {
                            if (local_ele > (T)0.) {
                                if (fp == 1) {
                                    denominator(ax, local_ind, local_ele, input, param.TOF, local_ele, TOFSum, DD, param.TOFCenters, param.sigma_x, D, param.nBins, lor, nRays, param.projType);
                                }
                                else if (fp == 2) {
                                    if (param.useMaskBP) {
                                        maskVal = param.maskBP[localIndX * d_N2 + localIndY * d_N3];
                                    }
                                    if (maskVal > 0)
                                        rhs(local_ele * temp, ax, local_ind, output, param.noSensImage, SensImage, local_ele, param.sigma_x, D, DD, param.TOFCenters, TOFSum, param.TOF, param.nBins, param.projType);
                                }
                            }
                        }
                        if (param.TOF)
                            D -= std::copysign(local_ele2, DD);
                        if (tempi < 0 || tempi >= param.Nx || tempj < 0 || tempj >= param.Ny || tempk < 0 || tempk >= param.Nz) {
                            break;
                        }
                        if (SPECT && (param.projType == 1)) {
                            L_SPECT += local_ele;
                        }
                    }
                    if (SPECT && (param.projType == 1)) {
                        temp /= L_SPECT;
                    }
                    if (param.projType > 1) {
                        if (ux < 0) {
                            for (int ii = tempi_a - 1; ii >= 0; ii--) {
                                int uu = orthDistance3D(ii, y_diff, x_diff, z_diff, center1[ii], center2, param.z_center, temp, tempj_a, tempk_a, xs, ys, detectors.zs, Nyx, kerroin, d_N1, d_N2, d_N3,
                                    param.Nz, param.bmin, param.bmax, param.Vmax, param.V, XY, ax, input, param.noSensImage, SensImage, output, local_ele, param.sigma_x, D, DD, param.TOFCenters, TOFSum, param.TOF, fp,
                                    param.projType, param.nBins, lor, nRays, param.useMaskBP, param.maskBP, attApu, SPECT, param.attenuationCorrection);
                                if (uu == 0)
                                    break;
                            }
                        }
                        else {
                            for (int ii = tempi_a + 1; ii < d_NNx; ii++) {
                                int uu = orthDistance3D(ii, y_diff, x_diff, z_diff, center1[ii], center2, param.z_center, temp, tempj_a, tempk_a, xs, ys, detectors.zs, Nyx, kerroin, d_N1, d_N2, d_N3,
                                    param.Nz, param.bmin, param.bmax, param.Vmax, param.V, XY, ax, input, param.noSensImage, SensImage, output, local_ele, param.sigma_x, D, DD, param.TOFCenters, TOFSum, param.TOF, fp,
                                    param.projType, param.nBins, lor, nRays, param.useMaskBP, param.maskBP, attApu, SPECT, param.attenuationCorrection);
                                if (uu == 0)
                                    break;
                            }
                        }
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