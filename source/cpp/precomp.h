#pragma once
#if defined(CUDA) || defined(HIP)
#include "cuda_error.h"
#define CTYPE3 float3
#define CTYPE2 float2
#elif defined(OPENCL)
#include "opencl_error.hpp"
#define CTYPE3 cl_float3
#define CTYPE2 cl_float2
#define ITYPE3 cl_int3
#elif defined(METAL)
#include "metal_error.hpp"
#include <simd/simd.h>
#define CTYPE3 simd::float3
#define CTYPE2 simd::float2
#define ITYPE3 simd::int3
#else
struct float3a {
	float x, y, z;
};
struct float2a {
	float x, y, z;
};
#define CTYPE3 float3a
#define CTYPE2 float2a
#endif
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <cstdint>
#include <fstream>
#include <limits>
#include <map>
#include <mutex>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#ifdef MATLAB
#include "mexFunktio.h"
#endif
#ifdef METAL
#define NS_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#include "Metal.hpp"
#endif

struct largeDimStruct {
	uint32_t NzOrig;
	float bzOrig;
	float bmaxZOrig;
	float d_Scale4ZOrig;
	int64_t imDimOrig;
	std::vector<uint32_t> Nz, NzPr, startPr, endPr;
	std::vector<float> bz;
	std::vector<float> bmaxZ;
	std::vector<float> d_Scale4Z;
	std::vector<int64_t> imDim, cumDim, imDimPr, cumDimPr;
};

// The dimensions of the input (src) and of the output (out) grid of the rotation of the SPECT
// rotation-based projector, together with the part of the output that is written and the part of the
// input that holds data, in dimensions 1 and 2. The two grids are the same volume without
// multi-resolution, but differ with it, since the rotated grid is always centered on the rotation axis
// while a subvolume is not. With every dimension set to the whole volume and both windows covering it
// the rotation behaves like a plain rotation of the entire image
// Dimension 0 (the depth) has a window of its own: outX and outNxW are the first index and the number of
// the planes of the output grid that are written (the output array is outNxW planes wide, the centering is
// on the whole grid, outNx), srcX and srcNxB the first index and the number of the planes of the source
// that hold data. Zero for outNxW / srcNxB means the whole of outNx / srcNx
typedef struct RotWindow_ {
	int srcNx = 0, srcNyFull = 0, srcNzFull = 0;
	int outNx = 0, outNyFull = 0;
	int outY = 0, outZ = 0, outNy = 0, outNz = 0;
	int srcY = 0, srcZ = 0, srcNy = 0, srcNz = 0;
	int outX = 0, outNxW = 0, srcX = 0, srcNxB = 0;
	// The dimensions of dimensions 0 and 1 that the output array is allocated with, i.e. its row and slice strides
	// (zero: outNxW and outNy, the array is as large as the window). The elements outside of the window are not written.
	// allocQ is set by the caller of rotateHelperType6 to have the output allocated with its three dimensions rounded up to a
	// multiple of it (zero: not), which makes the arrays of the views reuse the memory of the ones of the previous views
	int outSx = 0, outSy = 0, allocQ = 0;
} RotWindow;

// The scalars of the fused type-6 collimator-detector response kernels (cdrfAxial, cdrfTransFP,
// cdrfTransBP), see psfHelperType6 for their meaning. The response is worked out inside the kernels
typedef struct CDRFParams6_ {
	// Absolute plane index of the first plane of the arrays, the first visible plane and the radius of
	// the truncated Gaussians
	int xOff = 0, firstPlane = 0, radius = 0;
	// The weights of the response are made once per view and volume by cdrfWeights and read by the kernels (a large
	// region, or a radius that the kernels cannot make the weights for), instead of the kernels making them
	// in shared memory, see type6GlobalWeights
	bool globalW = false;
	// Distance from the first plane to the collimator surface (voxels) and the voxel sizes (mm)
	float dist = 0.f, dxv = 0.f, dyv = 0.f, dzv = 0.f;
	// Collimator resolution coefficients and intrinsic resolution, FWHM = sqrt((A * z + B)^2 + C^2)
	float Az = 0.f, Bz = 0.f, Axy = 0.f, Bxy = 0.f, C = 0.f;
	// FWHM to sigma and the normalization of the projection
	float sigmaCoeff = 0.f, norm = 1.f;
	// The variance (voxels^2 of dyv) that the rotation of the view adds transaxially (the interpolation of the rotation and, with
	// resampling onto the detector, of that), which the transaxial response is narrowed by (see rotVarType6, GridType6::rotVar). Negative:
	// the constant sigmaCoeff^2 that was used before the compensation depended on the view
	float comp = -1.f;
} CDRFParams6;

// The layout of the arrays of the axial kernel of the fused type-6 pipeline (cdrfAxial, cdrfAxialGW, cdrfAxialDirect). The arrays
// are logically N0 x N1 x N2 but may be allocated larger in the first two dimensions: inSx x inSy and outSx x outSy are the
// dimensions of the input and the output array (zero: N0 x N1). The attenuation is read from its own array, of the dimensions
// attSx x attSy x (anything) (zero: N0 x N1), at the offset (attX0, attY0, attZ0) of it. With attDx != 0 that array holds the
// accumulated attenuation (the sum along the depth, without the voxel size) and the attenuation factor is made of it as
// type6AttConv says (the average over the voxel by default)
typedef struct CDRFLayout6_ {
	int inSx = 0, inSy = 0, outSx = 0, outSy = 0;
	int attX0 = 0, attY0 = 0, attZ0 = 0, attSx = 0, attSy = 0;
	float attDx = 0.f;
} CDRFLayout6;

// How the fused type-6 kernels turn the accumulated attenuation S_k (the sum of the attenuation of the planes 0 ... k, including the plane
// k itself, along the depth from the detector side) into the attenuation factor of the voxels of plane k, voxel size dx. The environment
// variable TYPE6_ATTCONV: 2 (the default, variable not set or invalid) = exp(-dx * S_(k-1)) * (1 - exp(-t)) / t, t = dx * (S_k - S_(k-1)),
// the attenuation averaged exactly over the thickness of the emitting voxel (exact for a constant emission and attenuation within the
// voxel); it is the default because it is independent of the voxel size, which is required for a consistent multi-resolution (verified
// on CUDA: removes most of the multires attenuated-FP deficit); 0 = exp(-dx * S_k), the attenuation through the whole thickness of the
// voxel (the behaviour before 2026-10, ~1.5% lower attenuated projections, depends on the voxel size); 1 = exp(-dx * (S_(k-1) + S_k) / 2),
// the attenuation to the centre of the voxel (practically equal to 2). S_(-1) = 0. The ArrayFire pipeline (attenuationFactorType6) uses
// the same one. It is read once, when first used
inline int type6AttConv() {
	static const int mode = []() {
		int m = 2;
#ifdef _WIN32
		char* val = nullptr;
		size_t len = 0;
		if (_dupenv_s(&val, &len, "TYPE6_ATTCONV") == 0 && val != nullptr) {
			if (val[0] == '0' || val[0] == '1' || val[0] == '2')
				m = val[0] - '0';
			free(val);
		}
#else
		const char* val = std::getenv("TYPE6_ATTCONV");
		if (val != nullptr && (val[0] == '0' || val[0] == '1' || val[0] == '2'))
			m = val[0] - '0';
#endif
		return m;
	}();
	return mode;
}

// The per-view parameters of the batched fused type-6 kernels (rotateB, rotateBackB, cdrfAxialB, cdrfTransFPB, cdrfTransBPB,
// resampleToDetectorB, see TYPE6_BATCH in functions.hpp): one block of type6BatchStride floats per view of the batch, uploaded
// to the device once for every batch of views of a volume. The integers are stored exactly in the floats (all of them are below
// 2^24). The indices are those of the T6P_ defines of auxKernels.cl
constexpr int type6BatchStride = 53;
namespace type6P {
	enum {
		SRCNX = 0, SRCNY = 1, SRCNZ = 2, OUTNX = 3, OUTNYF = 4, COSA = 5, SINA = 6, OFFX = 7, OFFY = 8,
		OUTY = 9, OUTZ = 10, OUTNY = 11, OUTNZ = 12, SRCY = 13, SRCZ = 14, SRCNYB = 15, SRCNZB = 16,
		OUTX = 17, OUTNXW = 18, SRCX = 19, SRCNXB = 20,
		N0 = 21, N1 = 22, N2 = 23, ZSTART = 24, NZOUT = 25, XOFF = 26, FIRSTPLANE = 27, RADIUS = 28,
		ATTX0 = 29, ATTY0 = 30, ATTZ0 = 31, ATTDX = 32, DIST = 33,
		YSTART = 34, NOUT = 35, OUTSLOT = 36, VIEW = 37, PNY = 38, PNZ = 39, PY0 = 40, PZ0 = 41,
		TRA = 42, TRB = 43, AXA = 44, AXB = 45, NTR = 46, NAX = 47,
		Y0 = 48, NY = 49, IZ0 = 50, INZ = 51,
		COMP = 52
	};
}

// The transaxial variance that the bilinear interpolation of the rotation of a view adds to the image, compensated for by the
// transaxial response of the type-6 projector (cdrfSigma in auxKernels.cl: the response is narrowed by it). The rotation kernels (rotate,
// rotateB) give the element a of the rotated grid (centered coordinates, a = index - centre, row = transaxial) the value of the source at
// R a + off, R = [[cos, -sin], [sin, cos]], by the tent weights T(a, u) = tri(s_x - u_x) tri(s_y - u_y) of the source voxels u (centered on the
// source, u = index - centre) around s = R a + off. A source voxel u therefore spreads over the elements a around its own position
// q = R^T (u - off) of the grid, with a variance along the rows of V(u) = sum_a T (a_y - q_y)^2 / sum_a T. It depends on the angle and
// on where q is with respect to the lattice: with the lattices aligned (0, 90, 180, 270 degrees, integer offsets) it is zero, with
// fractional offsets f (1 - f), for any other angle it is 0 ... 1/4 with the mean 1/6 over the voxels (the second moment of the tent), and
// near the axis-aligned angles the lattices drift against each other across the volume. The response is made for the mean over the voxels
// of the volume (the activity is unknown), a single point source differs from it by 0.06 voxels^2 (standard deviation over the positions).
// The centre of the grid of n elements is that of the kernels, the exact centre n / 2 - 0.5 in the index coordinates (real division: the middle
// of the middle element for the odd n, the kernels compute (2 x + 1 - n) / 2 from it). The kernels (rotate, rotateB, rotateBackB), the depth
// window (depthWindowType6) and this function have to be changed together
inline double type6RotCentre(const int n) {
	return 0.5 * static_cast<double>(n) - 0.5;
}

// The mean of V(u) (see above) over the voxels of the srcNx x srcNy source, for the rotation cosa, sina (the values that the rotation kernel
// gets) with the offsets offX, offY (voxels) into the outNx x outNy rotated grid. Weighted by the total weight of each voxel (the mean of
// the second moment over all the contributions). A subset of the voxels is used (Halton points, rounded to voxels), all of them when it
// is smaller than the number of the samples. The weights do not depend on the axial direction. In double precision, a few thousand
// operations
inline double rotVarType6Compute(const double cosa, const double sina, const double offX, const double offY, const int srcNx, const int srcNy, const int outNx, const int outNy) {
	static constexpr int nSamples = 1024;
	// Exactly aligned, V is the same for every voxel
	const bool aligned = std::abs(sina) < 1e-6 || std::abs(cosa) < 1e-6;
	const double csx = type6RotCentre(srcNx), csy = type6RotCentre(srcNy);
	const double cox = type6RotCentre(outNx), coy = type6RotCentre(outNy);
	const int64_t nVox = static_cast<int64_t>(srcNx) * static_cast<int64_t>(srcNy);
	const bool all = nVox <= nSamples;
	const int nS = aligned ? 1 : (all ? static_cast<int>(nVox) : nSamples);
	// The Halton sequence (bases 2 and 3), made once
	static const std::vector<std::array<double, 2>> halton = []() {
		std::vector<std::array<double, 2>> h(nSamples);
		for (int i = 0; i < nSamples; i++) {
			double v[2];
			const int base[2] = { 2, 3 };
			for (int d = 0; d < 2; d++) {
				double f = 1.0, r = 0.0;
				int k = i + 1;
				while (k > 0) {
					f /= base[d];
					r += f * (k % base[d]);
					k /= base[d];
				}
				v[d] = r;
			}
			h[i] = { v[0], v[1] };
		}
		return h;
	}();
	const double hh = std::abs(cosa) + std::abs(sina);
	double num = 0., den = 0.;
	for (int ss = 0; ss < nS; ss++) {
		int ix, iy;
		if (aligned) {
			ix = srcNx / 2;
			iy = srcNy / 2;
		}
		else if (all) {
			ix = static_cast<int>(ss % srcNx);
			iy = static_cast<int>(ss / srcNx);
		}
		else {
			ix = (std::min)(static_cast<int>(halton[ss][0] * srcNx), srcNx - 1);
			iy = (std::min)(static_cast<int>(halton[ss][1] * srcNy), srcNy - 1);
		}
		const double ux = static_cast<double>(ix) - csx;
		const double uy = static_cast<double>(iy) - csy;
		// The position of the voxel in the rotated grid, a = R^T (u - off), and the elements that its tent weights reach
		const double qx = cosa * (ux - offX) + sina * (uy - offY);
		const double qy = -sina * (ux - offX) + cosa * (uy - offY);
		const int jx0 = static_cast<int>(std::ceil(qx + cox - hh)), jx1 = static_cast<int>(std::floor(qx + cox + hh));
		const int jy0 = static_cast<int>(std::ceil(qy + coy - hh)), jy1 = static_cast<int>(std::floor(qy + coy + hh));
		for (int jy = jy0; jy <= jy1; jy++) {
			const double ay = static_cast<double>(jy) - coy;
			for (int jx = jx0; jx <= jx1; jx++) {
				const double ax = static_cast<double>(jx) - cox;
				const double sx = cosa * ax - sina * ay + offX;
				const double sy = sina * ax + cosa * ay + offY;
				const double T = (std::max)(1. - std::abs(sx - ux), 0.) * (std::max)(1. - std::abs(sy - uy), 0.);
				num += T * (ay - qy) * (ay - qy);
				den += T;
			}
		}
	}
	return den > 0. ? num / den : 1. / 6.;
}

// The transaxial variance that the resampling of the projection of a volume onto the detector adds (the volumes that are not on the lattice
// of the detector, see GridType6::resample): the detector row r reads the rotated grid at the position A r + B by the tent weights, so a
// projected point at a lattice position i0 has the weights tri(A r + B - i0) over the rows, whose variance (voxels^2 of the grid) is
// sum_r w (A r + B - i0)^2 / sum_r w. Its mean over the lattice positions i0 of the nDet rows, weighted by the weights, which does not
// depend on the integer part of B. Zero for A = 1 and B an integer (the lattices coincide)
inline double detVarType6Compute(const int nDet, const double A, const double B) {
	double num = 0., den = 0.;
	for (int r = 0; r < nDet; r++) {
		const double t = A * static_cast<double>(r) + B;
		const double i0 = std::floor(t);
		for (int jj = 0; jj < 2; jj++) {
			const double d = t - (i0 + static_cast<double>(jj));
			const double w = (std::max)(1. - std::abs(d), 0.);
			num += w * d * d;
			den += w;
		}
	}
	return den > 0. ? num / den : 0.;
}

// rotVarType6Compute, kept for the rest of the run: it is asked for by every projection and backprojection of a view and volume, and the
// views (and the volumes) of the geometry are the same for all of them
inline double rotVarType6(const double cosa, const double sina, const double offX, const double offY, const int srcNx, const int srcNy, const int outNx, const int outNy) {
	static std::mutex mtx;
	static std::map<std::array<double, 8>, double> cache;
	const std::array<double, 8> key = { cosa, sina, offX, offY, static_cast<double>(srcNx), static_cast<double>(srcNy), static_cast<double>(outNx), static_cast<double>(outNy) };
	{
		std::lock_guard<std::mutex> lock(mtx);
		const auto it = cache.find(key);
		if (it != cache.end())
			return it->second;
	}
	const double v = rotVarType6Compute(cosa, sina, offX, offY, srcNx, srcNy, outNx, outNy);
	std::lock_guard<std::mutex> lock(mtx);
	cache[key] = v;
	return v;
}

// Whether the transaxial response is narrowed by the variance that the rotation of the view adds (the default), or by the constant
// sigmaCoeff^2 (= 0.18, a FWHM of one voxel) of before: the environment variable TYPE6_ROTCOMP = 0 gives the constant (for comparing).
// It is read once
inline bool type6RotComp() {
	static const bool on = []() {
		bool m = true;
#ifdef _WIN32
		char* val = nullptr;
		size_t len = 0;
		if (_dupenv_s(&val, &len, "TYPE6_ROTCOMP") == 0 && val != nullptr) {
			if (val[0] == '0')
				m = false;
			free(val);
		}
#else
		const char* val = std::getenv("TYPE6_ROTCOMP");
		if (val != nullptr && val[0] == '0')
			m = false;
#endif
		return m;
	}();
	return on;
}

// Test switch for the rotation images of the type-6 projector (the environment variable TYPE6_FORCE_TEX_FAIL), so that the paths of a
// device that cannot make them can be tested on one that can: 1 makes the probe of ProjectorClass::addProjector report that the images
// are not supported, and the projector falls back to buffers; 2 lets the probe pass but makes the creation of the first rotation image
// fail, which stops the projection with an error. Anything else (or not set) does nothing. It is read whenever a projector is created,
// not cached, so that it can be changed between runs in the same MATLAB/Python session
inline int type6ForceTexFailEnv() {
	int m = 0;
#ifdef _WIN32
	char* val = nullptr;
	size_t len = 0;
	if (_dupenv_s(&val, &len, "TYPE6_FORCE_TEX_FAIL") == 0 && val != nullptr) {
		if (val[0] == '1' || val[0] == '2')
			m = val[0] - '0';
		free(val);
	}
#else
	const char* val = std::getenv("TYPE6_FORCE_TEX_FAIL");
	if (val != nullptr && (val[0] == '1' || val[0] == '2'))
		m = val[0] - '0';
#endif
	return m;
}

typedef struct structForScalars {
	uint32_t projector_type = 1, attenuation_correction = 0, randoms_correction = 0, scatter = 0, normalization_correction = 0, 
		nColsD, nRowsD, nHeads = 1, size_z, subsets = 1, det_per_ring, Niter = 1, Nt = 1, subsetType = 0, nMultiVolumes = 0, nLayers = 1,
		nRekos = 1, osa_iter0 = 0, timestep0 = 0, nRekos2 = 0, subsetsUsed = 1, timestepsUsed = 1, TOFsubsets = 1, Nxy = 0U, NxOrig = 0U, NyOrig = 0U, NzOrig = 0U, NxPrior = 0U, NyPrior = 0U, NzPrior = 0U,
		BPType = 1, FPType = 1, adaptiveType = 0, rings = 0, FISTAType = 0, maskFPZ = 1, normZ = 1, maskBPZ = 1, currentSubset = 0;
	uint32_t platform = 0;
	std::vector<uint32_t> Nx{ 1, 0, 0, 0, 0, 0, 0 }, Ny{ 1, 0, 0, 0, 0, 0, 0 }, Nz{ 1, 0, 0, 0, 0, 0, 0 };
	float crystal_size_z = 0.f, epps = 1e-6f, sigma_x = 0.f, tube_width = 0.f, bmin = 0.f, bmax = 0.f, Vmax = 0.f, global_factor = 1.f,
		dL = 0.f, flat = 0.f, cylRadiusProj3 = 0.f, DSC = 0.f, helicalRadius = 0.f;//, T = 0.f
	std::vector<float> dx{ 0.f, 0.f }, dy{ 0.f, 0.f }, dz{ 0.f, 0.f }, bx{ 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f }, by{ 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f }, 
		bz{ 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f };
    bool useHalf = false;
	bool use_psf = false, TOF = false, SPECT = false, pitch = false, PET = false, meanFP = false, meanBP = false,
		maskFP = false, maskBP = false, orthXY = false, orthZ = false, CT = false, atomic_64bit = false, atomic_32bit = false, loadTOF = true,
		saveIter = false, enforcePositivity = false, computeSensImag = false, useMAD = true, useImages = false, eFOV = false,
		useExtendedFOV = false, use64BitIndices = false, TGV2D = false, multiResolution = false, offset = false, offsetSqrt = false, relaxScaling = false,
		computeRelaxation = false, storeFP = false, deconvolution = false, CTAttenuation = true, largeDim = false, usingLinearizedData = false, storeResidual = false,
		useBuffers = true, useFDKWeights = false, indexBased = false, FISTAAcceleration = false, stochastic = false, useTotLength = true,
		useParallelBeam = false, useHelical = false;
	// Set this to false to explicitly disable fastPDHG
	// On MATLAB/Octave, you can do this just with options.fastPDHG = false;
	bool fastPDHG = true;
	// True only when Nt > 1 AND some sensitivity-affecting input (image-based attenuation or the
	// multiplicative corrVector/extraCorr) genuinely differs per timestep. Computed once, early, in
	// reconstructionAF() and consumed there plus in functions.hpp (transferControl, backprojectionType6)
	// and subiterStep.h (computeOSEstimates) to pick vec.Summ's timestep slot: sIdx = sensPerTimestep ?
	// timestep : 0. When false (the common case, including all Nt == 1 runs), every timestep shares the
	// single sensitivity image computed from timestep 0 -- see reconstructionAF.h for why per-timestep
	// slots are otherwise unnecessary (normalization and measurement-domain attenuation are single
	// buffers shared by all timesteps).
	bool sensPerTimestep = false;
	int64_t Nf = 0;
	std::vector<CTYPE3> d_Scale;
	std::vector<CTYPE3> d_Scale4;
	std::vector<CTYPE2> dSize;
	CTYPE2 dSizeBP;
	uint8_t raw = 0, fp = 0, listmode = 0;
	int8_t verbose = 0;
	uint16_t n_rays = 1, n_rays3D = 1;
	uint32_t g_dim_x = 0u, g_dim_y = 0u, g_dim_z = 0u, deblur_iterations = 0U;
	// User-defined local/work-group (block) size per dimension. Any negative value means "use the
	// built-in default" for that dimension (see ProjectorClass::addProjector).
	int32_t localSize[3] = { -1, -1, -1 };
	// Compute the (spatial) regularization/prior only every regEveryIter-th time it is reached. A value
	// of 1 (or less) computes it every time (original behavior). The first and last time are always
	// computed. See subiterStep.h and iterStep.h.
	int32_t regEveryIter = 1;
	int64_t nBins = 1, nProjections = 0, numelY = 0, numelZ = 0, TOFSize = 0, seed = -1;
	std::vector<int64_t> im_dim{ 1 };
	size_t size_of_x, size_atten = 1, size_norm = 1, size_center_x, size_center_y, size_center_z, size_V = 1, size_scat = 1, kokoTOF = 0, kokoNonTOF = 0, sizeLOR,
		sizeL, sizeXY, sizeZ, saveIterationsMiddle = 0ULL;
	// Variables to check for correct size
	size_t size_maskFP = 0ULL, size_maskBP = 0ULL, size_meas = 0ULL;
	uint32_t* saveNIter = nullptr;
	float* T = nullptr;
	float* V = nullptr;
	float* x_center = nullptr, *y_center = nullptr, *z_center = nullptr;
	float* gaussian = nullptr;
	float* TOFCenter = nullptr;
	uint64_t* pituus, length;
	std::vector<uint32_t> usedDevices;
	largeDimStruct lDimStruct;
	float coneOfResponseStdCoeffA = 0.01f, coneOfResponseStdCoeffB = 0.01f, coneOfResponseStdCoeffC = 0.01f;
	float ellipseCenterX = 0.f, ellipseCenterY = 0.f, ellipseCenterZ = 0.f, ellipseRadiusX = 1.f, ellipseRadiusY = 1.f, ellipseRadiusZ = 1.f,
		ellipsePower = std::numeric_limits<float>::max(); // Box support, see extendRayToEllipse
	// Compute the collimator-detector response of the rotation-based SPECT projector separately for
	// every projection instead of shifting a precomputed filter bank
	bool CDRFOnTheFly = false;
	// The model of the collimator response of the on-the-fly CDRF (option CDRFModel): 0 = Gaussian (default), 1 = square holes (a triangle
	// convolved with the intrinsic Gaussian), see cdrfTap in auxKernels.cl
	uint32_t cdrfModel6 = 0U;
	// How the backprojection rotates a view back into the image (option type6BPRotation): 0 = the exact adjoint of the bilinear rotation of the
	// forward projection (default), 1 = the fast bilinear gather at the position of the inverse rotation (-DTYPE6_BP_FAST, see rotate in auxKernels.cl)
	uint32_t bpRot6 = 0U;
	// Depth-dependent collimator-detector response coefficients for projector type 6: {Az, Bz, Axy, Bxy, C}
	// with FWHM_axis(z_mm) = sqrt((A_axis * z + B_axis)^2 + C^2), z the distance from the collimator surface in mm
	float CDRFCoeffs6[5] = { 0.f, 0.f, 0.f, 0.f, 0.f };
	// The fraction of the available device memory that the type-6 projector may use to keep the attenuation of the
	// forward projection of a subset for its backprojection (option type6AttenuationCache, clamped to [0, 0.9], 0 = keep nothing)
	float type6AttenuationCache = 0.4f;
	// The largest absolute lateral shift of the detector panel over all projections of the rotation-based
	// SPECT projector, in mm (the maximum of |blurPlanes2F| * dx[0]). Used to size the rotated grid of the
	// main volume with multi-resolution, see rotGridType6
	float maxPanelShift6 = 0.f;
	// The largest number of slices that the rotated grid of every volume is padded with at either axial end with
	// multi-resolution (the radius bound), see computePadZType6. Empty (i.e. no padding) without multi-resolution
	std::vector<int> padZ6;
	// The number of slices that the rotated grid of every volume is actually padded with at its low and high axial end,
	// at most padZ6: only as many as the detector extends past that edge of the volume, see computePadZType6
	std::vector<int> padZLo6, padZHi6;
} scalarStruct;

// The position of the rotation axis of the gantry (x, y) and of the axial center of the detector (z) in
// the coordinates of the volumes. The volumes are placed with respect to the absolute origin (the center
// of the main volume is at the object offset), and so is the SPECT detector geometry: the radius of
// rotation and the distance from the axis to the detector are relative to the rotation axis, which is the
// origin, and the detector columns are centered on z = 0 (SPECTDetectorCoordinates.m). The center of the main
// volume is therefore the axis and the detector center only when there is no object offset
inline void axisType6(const scalarStruct&, float& ax, float& ay, float& az) {
	ax = 0.f;
	ay = 0.f;
	az = 0.f;
}

// Whether the main volume is centered on the rotation axis transaxially (to within 1e-3 of a voxel),
// i.e. whether it can be rotated into itself
inline bool mainCenteredXYType6(const scalarStruct& inputScalars) {
	float ax = 0.f, ay = 0.f, az = 0.f;
	axisType6(inputScalars, ax, ay, az);
	const float cx = inputScalars.bx[0] + static_cast<float>(inputScalars.Nx[0]) * inputScalars.dx[0] * .5f;
	const float cy = inputScalars.by[0] + static_cast<float>(inputScalars.Ny[0]) * inputScalars.dy[0] * .5f;
	return std::abs(cx - ax) < 1e-3f * inputScalars.dx[0] && std::abs(cy - ay) < 1e-3f * inputScalars.dy[0];
}

// Whether the main volume is centered on the rotation axis and on the axial center of the detector
inline bool mainCenteredType6(const scalarStruct& inputScalars) {
	float ax = 0.f, ay = 0.f, az = 0.f;
	axisType6(inputScalars, ax, ay, az);
	const float cz = inputScalars.bz[0] + static_cast<float>(inputScalars.Nz[0]) * inputScalars.dz[0] * .5f;
	return mainCenteredXYType6(inputScalars) && std::abs(cz - az) < 1e-3f * inputScalars.dz[0];
}

// Whether the main volume is rotated into an axis-centered disk grid instead of into itself, which is the
// case with multi-resolution and with a main volume that is off-centre transaxially
inline bool usesDiskGridType6(const scalarStruct& inputScalars) {
	return inputScalars.nMultiVolumes > 0 || !mainCenteredXYType6(inputScalars);
}

// Whether the rotated grids are padded axially: with multi-resolution and with a main volume that is off-centre
// in any direction. The detector is centered on z = 0, so with an axial offset its center is not the volume
// center and the axial blur past the edge of the volume would be truncated. The grid itself is still the
// volume when the main volume is centered transaxially, see usesDiskGridType6
inline bool usesPadZType6(const scalarStruct& inputScalars) {
	return inputScalars.nMultiVolumes > 0 || !mainCenteredType6(inputScalars);
}

// The axial padding bound of the rotated grid of volume ii, in slices of that volume: the largest padding that either
// end of the grid can have (the axial radius of the response). Zero without multi-resolution
inline int padZType6(const scalarStruct& inputScalars, const int ii) {
	if (!usesPadZType6(inputScalars) || ii < 0 || static_cast<size_t>(ii) >= inputScalars.padZ6.size())
		return 0;
	return inputScalars.padZ6[ii];
}

// The axial padding of the rotated grid of volume ii at its low (lo) and high (hi) end, in slices of that volume. The grid of
// volume ii is lo + Nz[ii] + hi slices deep and its slice k is slice k - lo of the volume. Both are zero without
// multi-resolution, in which case the grid is the volume itself
inline void padZLoHiType6(const scalarStruct& inputScalars, const int ii, int& lo, int& hi) {
	lo = hi = 0;
	if (!usesPadZType6(inputScalars) || ii < 0 || static_cast<size_t>(ii) >= inputScalars.padZLo6.size() || static_cast<size_t>(ii) >= inputScalars.padZHi6.size())
		return;
	lo = inputScalars.padZLo6[ii];
	hi = inputScalars.padZHi6[ii];
}
inline int padZLoType6(const scalarStruct& inputScalars, const int ii) {
	int lo = 0, hi = 0;
	padZLoHiType6(inputScalars, ii, lo, hi);
	return lo;
}
inline int padZHiType6(const scalarStruct& inputScalars, const int ii) {
	int lo = 0, hi = 0;
	padZLoHiType6(inputScalars, ii, lo, hi);
	return hi;
}

// The rotation-based SPECT projector rotates every volume around the rotation axis of the gantry,
// which is the origin (see axisType6). Only a main volume that is centered on it can be rotated into itself, so a
// multi-resolution subvolume (or an off-centre main volume) is rotated into a grid that is centered on the axis instead and is large
// enough to hold the subvolume at any angle, that is the disk swept by its corners. The grid keeps the
// voxel size of the subvolume, so the rotation itself never changes the sampling.
// Without multi-resolution and with the main volume centered on the axis, the main volume is rotated into
// itself, exactly as it always has been. With multi-resolution or an off-centre main volume the main volume is small and holds body activity close to its corners, which the
// rotation and the panel shift (that moves the panel center to the center of the grid) would move out of
// a grid of the same size, so it uses the same axis-centered disk grid as the subvolumes (R is then half
// the diagonal of the main volume). When the detector plane is wider than the disk (an uncropped
// projection padded to the extended FOV), the plane is placed on the grid with the panel shift, so the
// transaxial dimension is widened by the largest panel shift (maxPanelShift6) to keep the whole shifted
// disk inside the grid. The depth dimension is not affected by the lateral shift.
// The main volume is centered on the axis, so its disk grid follows the parity (even/odd) of the main
// volume dimensions Nx0 and Ny0, which keeps the detector window of the cropped projections (whose size
// has the parity of Ny0) on an integer offset in the grid. The subvolumes always use even dimensions
inline void rotGridType6(const scalarStruct& inputScalars, const int ii, int& outNx, int& outNy) {
	if (ii == 0 && inputScalars.nMultiVolumes == 0 && mainCenteredXYType6(inputScalars)) {
		outNx = static_cast<int>(inputScalars.Nx[0]);
		outNy = static_cast<int>(inputScalars.Ny[0]);
		return;
	}
	float cx = 0.f, cy = 0.f, cz = 0.f;
	axisType6(inputScalars, cx, cy, cz);
	const float x0 = inputScalars.bx[ii] - cx;
	const float x1 = x0 + static_cast<float>(inputScalars.Nx[ii]) * inputScalars.dx[ii];
	const float y0 = inputScalars.by[ii] - cy;
	const float y1 = y0 + static_cast<float>(inputScalars.Ny[ii]) * inputScalars.dy[ii];
	const float R = std::sqrt((std::max)(x0 * x0, x1 * x1) + (std::max)(y0 * y0, y1 * y1));
	// Smallest dimension that covers 2 * r at the voxel size d, with the parity of the main volume for the
	// main volume and even otherwise. Even dimensions keep the center of the grid exactly on the rotation
	// axis, and so does the parity of the main volume, which is centered on it
	const auto fit = [ii](const float r, const float d, const uint32_t Nmain) {
		int n = static_cast<int>(std::ceil(2.f * r / d));
		const int par = ii == 0 ? static_cast<int>(Nmain % 2U) : 0;
		if ((n % 2) != par)
			n++;
		return n;
	};
	float Rny = R;
	// A detector plane that is wider than the disk. The projections are on the lattice of the main volume
	if (0.5f * static_cast<float>(inputScalars.nRowsD) * inputScalars.dy[0] > R)
		Rny = R + inputScalars.maxPanelShift6;
	outNx = fit(R, inputScalars.dx[ii], inputScalars.Nx[0]);
	outNy = fit(Rny, inputScalars.dy[ii], inputScalars.Ny[0]);
}

// How far (in voxels, not yet rounded up) the response of a plane reaches, where it has fallen to 1e-6 of its peak (see cdrfRadius in
// auxKernels.cl): Rg is the geometric FWHM and C the intrinsic one (mm), dd the voxel size of the direction (mm). The Gaussian model
// (0) is 5.2565 sigma of the Gaussian of FWHM sqrt(Rg^2 + C^2), the square model (1) the half base of the triangle (Rg) plus 5.2565 sigma
// of the intrinsic Gaussian. The compensation of the blur of the rotation only narrows the response, so it is left out
inline float cdrfExtentType6(const uint32_t model, const float Rg, const float C, const float dd, const float sigmaCoeff) {
	if (model == 1U)
		return Rg / dd + 5.2565f * (C / dd * sigmaCoeff);
	return 5.2565f * (std::sqrt(Rg * Rg + C * C) / dd * sigmaCoeff);
}

// With multi-resolution (or an off-centre main volume) the volumes do not span the whole axial FOV, and the axial blur of the
// collimator-detector response of activity close to the axial edge of a volume spills over that edge onto
// detector columns that are beyond it. The response is applied on the rotated grid, so the grid is padded
// axially by up to padZ6[ii] slices at both ends (grid slice k is slice k - padZLo6[ii] of the volume) to hold that
// spill. The padding is only needed at an axial edge that lies inside the axial coverage of the detector (the nColsD columns
// at the voxel pitch of the main volume, centered on z = az of axisType6): beyond the detector there is no column to spill onto.
// So the padding at an end is the number of slices outside the volume that the tent interpolation of the detector resampling
// reaches from the outermost detector column centre (half a main-volume pitch inside the detector edge), at most padZ6, at least zero. Otherwise the blur is truncated at the edge of the volume, and the counts and the centroid near
// the edge are wrong.
// The padding is a fixed value per volume, so that the array shapes stay the same for every projection.
// It is the largest axial radius of the response that the projector can use for the volume, that is the
// one of psfHelperType6 (5.2565 sigma, rounded up to a multiple of radiusStep) with the axial sigma at the
// largest depth of any plane of the grid over all projections. The response is monotonic in the depth,
// so the widest one is at either end of the range and the depth of the collimator surface is the other
// end. The depth of the last plane of the grid is (grid.Nx - 0.5 - dist) * dx[ii] with the distance dist
// of gridType6, which is the smallest for the smallest distFloat (the distance of the collimator surface
// in voxels of the main volume, see gridType6). Zero without multi-resolution
inline void computePadZType6(scalarStruct& inputScalars, const float* distFloat, const int64_t nViews = 0) {
	inputScalars.padZ6.assign(static_cast<size_t>(inputScalars.nMultiVolumes) + 1U, 0);
	inputScalars.padZLo6.assign(static_cast<size_t>(inputScalars.nMultiVolumes) + 1U, 0);
	inputScalars.padZHi6.assign(static_cast<size_t>(inputScalars.nMultiVolumes) + 1U, 0);
	// The geometry arrays hold nViews entries (all the frames of dynamic data), or one for each of the projections without it
	const int64_t nViewsUsed = nViews > 0 ? nViews : inputScalars.nProjections;
	if (!usesPadZType6(inputScalars) || distFloat == nullptr || nViewsUsed <= 0)
		return;
	float minDist = distFloat[0];
	for (int64_t pp = 1; pp < nViewsUsed; pp++)
		minDist = (std::min)(minDist, distFloat[pp]);
	const float Az = inputScalars.CDRFCoeffs6[0], Bz = inputScalars.CDRFCoeffs6[1];
	const float C = inputScalars.CDRFCoeffs6[4];
	const float sigmaCoeff = 1.f / (2.f * std::sqrt(2.f * std::log(2.f)));
	const int radiusStep = 8;
	const float collimator = (minDist - .5f * static_cast<float>(inputScalars.Nx[0])) * inputScalars.dx[0];
	for (int ii = 0; ii <= static_cast<int>(inputScalars.nMultiVolumes); ii++) {
		int gridNx = 0, gridNy = 0;
		rotGridType6(inputScalars, ii, gridNx, gridNy);
		const float dist = (collimator + .5f * static_cast<float>(gridNx) * inputScalars.dx[ii]) / inputScalars.dx[ii];
		const float zMax = (static_cast<float>(gridNx) - .5f - dist) * inputScalars.dx[ii];
		float extentMax = 0.f;
		for (int ee = 0; ee < 2; ee++) {
			const float zz = ee == 0 ? 0.f : (std::max)(zMax, 0.f);
			const float Rg = (std::max)(Az * zz + Bz, 0.f);
			extentMax = (std::max)(extentMax, cdrfExtentType6(inputScalars.cdrfModel6, Rg, C, inputScalars.dz[ii], sigmaCoeff));
		}
		const int r = static_cast<int>(std::ceil(extentMax));
		const int padZ = ((r + radiusStep - 1) / radiusStep) * radiusStep;
		inputScalars.padZ6[ii] = padZ;
		// The tent interpolation of the detector resampling reaches, for the outermost detector column centres (half a pitch of
		// the main volume inside the detector edges, which centered on az spans nColsD columns at the pitch dz[0]), the slices
		// of the volume that surround them. A column centre that is outside the first (last) slice centre of this volume, which
		// happens on a coarser lattice even when the edges coincide, reaches the slice outside it, which holds the blur spill. So
		// the padding is the number of slices outside the volume that the support reaches, at most padZ. A column centre that
		// is on a slice centre exactly (volumes on the detector lattice) needs none
		float ax = 0.f, ay = 0.f, az = 0.f;
		axisType6(inputScalars, ax, ay, az);
		const float zColLo = az - .5f * static_cast<float>(inputScalars.nColsD) * inputScalars.dz[0] + .5f * inputScalars.dz[0];
		const float zColHi = az + .5f * static_cast<float>(inputScalars.nColsD) * inputScalars.dz[0] - .5f * inputScalars.dz[0];
		const float posLo = (zColLo - inputScalars.bz[ii]) / inputScalars.dz[ii] - .5f;
		const float posHi = (zColHi - inputScalars.bz[ii]) / inputScalars.dz[ii] - .5f;
		const float eps = 1e-4f;
		const int needLo = -static_cast<int>(std::floor(posLo + eps));
		const int needHi = static_cast<int>(std::ceil(posHi - eps)) - (static_cast<int>(inputScalars.Nz[ii]) - 1);
		inputScalars.padZLo6[ii] = (std::max)((std::min)(padZ, needLo), 0);
		inputScalars.padZHi6[ii] = (std::max)((std::min)(padZ, needHi), 0);
	}
}

// The attenuation-map source lattice used for volume ii of the rotation-based SPECT projector. The map
// is read in the rotated frame with the same y/z lattice as the emission grid (grid.Ny, grid.Nz), but
// extended along x, toward the detector, so that the path through every volume's material is
// accumulated even though each volume is rotated and attenuated separately.
//
// With a single volume (nMultiVolumes == 0) the map is the main volume as is: attNx == grid.Nx and the
// source lattice is (Nx0, Ny0, Nz0), i.e. identical to the current (pre-multires) behavior.
//
// With multiple volumes the attenuation source has to be large enough to hold the union bounding box of
// every volume, at volume ii's own voxel size, centered on the rotation axis (the center of the main
// volume). Runion is the largest distance from the axis to a corner of that union box, in the xy plane.
// attNx/srcNy are the smallest dimensions covering 2*Runion at volume ii's spacing, adjusted up by at
// most one voxel to keep the same parity (even/odd) as the emission grid, so the two lattices share a
// common center voxel/boundary
inline void attGridType6(const scalarStruct& inputScalars, const int ii, const int gridNx, const int gridNy, int& attNx, int& srcNx, int& srcNy) {
	// A single main volume that is centered on the axis is used as it is. An off-centre one is rotated into a
	// disk grid, and its map is resampled onto an axis-centered lattice like that of the multi-resolution volumes
	if (inputScalars.nMultiVolumes == 0 && mainCenteredXYType6(inputScalars)) {
		attNx = gridNx;
		srcNx = gridNx;
		srcNy = gridNy;
		return;
	}
	float cx = 0.f, cy = 0.f, cz = 0.f;
	axisType6(inputScalars, cx, cy, cz);
	float Runion = 0.f;
	for (int jj = 0; jj <= static_cast<int>(inputScalars.nMultiVolumes); jj++) {
		const float x0 = inputScalars.bx[jj] - cx;
		const float x1 = x0 + static_cast<float>(inputScalars.Nx[jj]) * inputScalars.dx[jj];
		const float y0 = inputScalars.by[jj] - cy;
		const float y1 = y0 + static_cast<float>(inputScalars.Ny[jj]) * inputScalars.dy[jj];
		const float R = std::sqrt((std::max)(x0 * x0, x1 * x1) + (std::max)(y0 * y0, y1 * y1));
		Runion = (std::max)(Runion, R);
	}
	const int minNx = static_cast<int>(std::ceil(2.f * Runion / inputScalars.dx[ii]));
	const int minNy = static_cast<int>(std::ceil(2.f * Runion / inputScalars.dy[ii]));
	// Bump the size up by one voxel, if needed, to match the parity of the emission grid
	attNx = (((minNx - gridNx) % 2) != 0) ? minNx + 1 : minNx;
	srcNx = attNx;
	srcNy = (((minNy - gridNy) % 2) != 0) ? minNy + 1 : minNy;
}

// The largest rotation grid over every volume, which is how large the images used by the rotation have
// to be. Both the volume itself (forward projection) and its rotated grid (backprojection) are
// uploaded into them, as is the attenuation source lattice (attGridType6) of every volume
inline void rotGridMaxType6(const scalarStruct& inputScalars, int& maxNx, int& maxNy, int& maxNz) {
	maxNx = maxNy = maxNz = 0;
	for (int ii = 0; ii <= static_cast<int>(inputScalars.nMultiVolumes); ii++) {
		int outNx = 0, outNy = 0;
		rotGridType6(inputScalars, ii, outNx, outNy);
		maxNx = (std::max)(maxNx, (std::max)(outNx, static_cast<int>(inputScalars.Nx[ii])));
		maxNy = (std::max)(maxNy, (std::max)(outNy, static_cast<int>(inputScalars.Ny[ii])));
		// The volume is uploaded at the offset of the axial padding of its grid, and the grid is that
		// much deeper than the volume
		maxNz = (std::max)(maxNz, static_cast<int>(inputScalars.Nz[ii]) + padZLoType6(inputScalars, ii) + padZHiType6(inputScalars, ii));
		int attNx = 0, srcNx = 0, srcNy = 0;
		attGridType6(inputScalars, ii, outNx, outNy, attNx, srcNx, srcNy);
		maxNx = (std::max)(maxNx, (std::max)(attNx, srcNx));
		maxNy = (std::max)(maxNy, srcNy);
	}
}

// The rotation-source images of the type-6 projector (ProjectorClass::d_rotI). Slots 0 and 1 are the ones of the
// volume-by-volume projector (0 = emission volume or the region that is rotated back, 1 = attenuation source), both as large
// as rotGridMaxType6. The projector that handles all volumes of a view together (forwardProjectionType6All,
// backprojectionType6All) keeps the emission image of every volume (slot rotSlotEmissionType6 + volume) and the attenuation
// source of every attenuation group (slot rotSlotAttenuationType6 + group, see attGroupsType6) in a slot of its own, each
// created, lazily, with the size of the image that it holds. There are at most 7 volumes and so at most 7 groups
constexpr int rotSlotEmissionType6 = 2;
constexpr int rotSlotAttenuationType6 = 2 + 7;
constexpr int rotSlotsType6 = 2 + 7 + 7;

// Volumes that have the same attenuation source lattice share one attenuation rotation and accumulation per view. The
// lattice is that of attenuationSourceLatticeType6: the axis-centered (srcNx x srcNy) lattice at the voxel size of the volume,
// with attNx planes, and its z lattice. Two volumes are in the same group when their voxel sizes and their (attNx, srcNx,
// srcNy) are the same and their z lattices are aligned (the offset of their first slices is a whole number of voxels), in
// which case the lattice of the group is the union of theirs: the source is the map resampled on it once, and the
// volume's own part of it is what the volume would have got on its own lattice (the axial padding of a rotated grid,
// padZLoHiType6, is zero in the source, and nothing reads the attenuation there, since the emission volume has no data in it
// and the backprojection does not keep it). The rotated attenuation is the one of the group frame: the rows of a member's
// grid (grid Ny rows, centered on the axis) are the rows (G - gridNy) / 2 + r of the group, G the largest grid Ny of the
// members, and the planes of a member's grid (grid Nx planes) the planes (attNx - gridNx) / 2 + k. The parity of the
// grids of the members is the same as that of the group (attGridType6 makes attNx and srcNy match the grid), which is what
// makes the centers agree; it is checked anyway
typedef struct AttGroupType6_ {
	std::vector<int> members;
	// The lattice of the source
	int attNx = 0, srcNx = 0, srcNy = 0;
	// The largest transaxial dimension of the rotated grid of the members
	int rows = 0;
	// The axial lattice of the source: the first slice is at bz, nz slices; pad is the largest axial padding of the grids of the
	// members (the source is placed padLo slices into the group frame, which is padLo + nz + padHi slices deep, padLo and padHi being the
	// largest padding at the low and at the high end of the members), and zOff[ii] the slice
	// of the lattice that the first slice of volume ii is
	float bz = 0.f;
	int nz = 0, padLo = 0, padHi = 0;
	std::vector<int> zOff;
} AttGroupType6;

inline std::vector<AttGroupType6> attGroupsType6(const scalarStruct& inputScalars) {
	std::vector<AttGroupType6> groups;
	const int nV = static_cast<int>(inputScalars.nMultiVolumes) + 1;
	for (int ii = 0; ii < nV; ii++) {
		int gNx = 0, gNy = 0;
		rotGridType6(inputScalars, ii, gNx, gNy);
		int attNx = 0, srcNx = 0, srcNy = 0;
		attGridType6(inputScalars, ii, gNx, gNy, attNx, srcNx, srcNy);
		int found = -1;
		for (size_t gg = 0; gg < groups.size() && found < 0; gg++) {
			const int jj = groups[gg].members[0];
			if (inputScalars.dx[ii] != inputScalars.dx[jj] || inputScalars.dy[ii] != inputScalars.dy[jj] || inputScalars.dz[ii] != inputScalars.dz[jj]
				|| attNx != groups[gg].attNx || srcNx != groups[gg].srcNx || srcNy != groups[gg].srcNy)
				continue;
			int jNx = 0, jNy = 0;
			rotGridType6(inputScalars, jj, jNx, jNy);
			if (((gNx - jNx) % 2) != 0 || ((gNy - jNy) % 2) != 0)
				continue;
			const float off = (inputScalars.bz[ii] - inputScalars.bz[jj]) / inputScalars.dz[ii];
			if (std::abs(off - std::round(off)) < 1e-3f)
				found = static_cast<int>(gg);
		}
		if (found < 0) {
			AttGroupType6 grp;
			grp.attNx = attNx;
			grp.srcNx = srcNx;
			grp.srcNy = srcNy;
			grp.zOff.assign(static_cast<size_t>(nV), 0);
			groups.push_back(grp);
			found = static_cast<int>(groups.size()) - 1;
		}
		groups[found].members.push_back(ii);
	}
	for (size_t gg = 0; gg < groups.size(); gg++) {
		AttGroupType6& grp = groups[gg];
		float bz = inputScalars.bz[grp.members[0]];
		for (size_t mm = 1; mm < grp.members.size(); mm++)
			bz = (std::min)(bz, inputScalars.bz[grp.members[mm]]);
		grp.bz = bz;
		for (size_t mm = 0; mm < grp.members.size(); mm++) {
			const int ii = grp.members[mm];
			int gNx = 0, gNy = 0;
			rotGridType6(inputScalars, ii, gNx, gNy);
			grp.rows = (std::max)(grp.rows, gNy);
			grp.padLo = (std::max)(grp.padLo, padZLoType6(inputScalars, ii));
			grp.padHi = (std::max)(grp.padHi, padZHiType6(inputScalars, ii));
			grp.zOff[ii] = static_cast<int>(std::lround((inputScalars.bz[ii] - bz) / inputScalars.dz[ii]));
			grp.nz = (std::max)(grp.nz, grp.zOff[ii] + static_cast<int>(inputScalars.Nz[ii]));
		}
	}
	return groups;
}

#ifdef OPENCL
// Struct for the various estimates used in OpenCL kernels
typedef struct _OpenCL_im_vectors {
	cl::Buffer d_meanFP, d_meanBP;
	cl::Buffer	d_im;
	std::vector<cl::Buffer> d_rhs_os;
	cl::Image3D d_image_os, d_image_os_int;
} OpenCL_im_vectors;
#elif defined(CUDA) || defined(HIP)
typedef struct _CUDA_im_vectors {
	CUdeviceptr d_meanFP, d_meanBP;
	CUdeviceptr* d_im;
	std::vector<CUdeviceptr*> d_rhs_os;
	CUtexObject d_image_os, d_image_os_int;
} CUDA_im_vectors;
#elif defined(METAL)
typedef struct _METAL_im_vectors {
	NS::SharedPtr<MTL::Buffer> d_im;
    std::vector<NS::SharedPtr<MTL::Buffer>> d_rhs_os;
    NS::SharedPtr<MTL::Texture> d_image_os, d_image_os_int;
} METAL_im_vectors;
#else
struct CPUVectors {
	float* d_meanFP, *d_meanBP, *d_im_os;
	std::vector<float*> d_rhs_os;
};
#endif

// Optional device-resident inputs/outputs for implementation 5 (e.g. MATLAB gpuArray data)
// Every member is a device pointer valid in the backend's current context; nullptr means the
// corresponding host array is used instead, i.e. the original behavior
typedef struct _deviceIO {
	const void* im = nullptr;      // Forward projection input image (type 1)
	const void* meas = nullptr;    // Backprojection input measurements (type 2)
	void* output = nullptr;        // FP measurement output (type 1) or BP image output (type 2)
	void* sensIm = nullptr;        // BP sensitivity image output (type 2, only when no_norm == 0)
} deviceIO;

inline void mexPrint(const char* str) {
#ifdef MATLAB
	mexPrintf("%s\n", str);
	mexEvalString("pause(.0001);");
#else
	fprintf(stdout, "%s\n", str);
	fflush(stdout);
#endif
}

inline void mexPrintVar(const char* str, const int var) {
#ifdef MATLAB
	mexPrintf("%s%d\n", str, var);
	mexEvalString("pause(.0001);");
#else
	fprintf(stdout, "%s%d\n", str, var);
	fflush(stdout);
#endif
}

inline void mexPrintVarf(const char* str, const float var) {
#ifdef MATLAB
	mexPrintf("%s%f\n", str, var);
	mexEvalString("pause(.0001);");
#else
	fprintf(stdout, "%s%f\n", str, var);
	fflush(stdout);
#endif
}

template <typename T>
inline void mexPrintBase(const char* str, const T var) {
#ifdef MATLAB
	mexPrintf(str, var);
#else
	fprintf(stdout, str, var);
#endif
}

template <typename T, typename K>
inline void mexPrintBase(const char* str, const T var1, const K var2) {
#ifdef MATLAB
	mexPrintf(str, var1, var2);
#else
	fprintf(stdout, str, var1, var2);
#endif
}

template <typename T, typename K, typename C>
inline void mexPrintBase(const char* str, const T var1, const K var2, const C var3) {
#ifdef MATLAB
	mexPrintf(str, var1, var2, var3);
#else
	fprintf(stdout, str, var1, var2, var3);
#endif
}

inline void mexEval() {
#ifdef MATLAB
	mexEvalString("pause(.0001);");
#else
	fflush(stdout);
#endif
}

inline void mexWarning(const char* str) {
#ifdef MATLAB
	mexWarnMsgTxt(str);
#endif
}
