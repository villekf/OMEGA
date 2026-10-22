/**************************************************************************
* Header for ArrayFire functions. Used by implementation 2.
*
* Copyright(C) 2019-2024 Ville - Veikko Wettenhovi
*
* This program is free software : you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <https://www.gnu.org/licenses/>.
***************************************************************************/
#pragma once
#define NOMINMAX
#if defined(CUDA) || defined(HIP)
#include "ProjectorClass.h"
#if defined(HIP)
// hipDeviceptr_t (== CUdeviceptr here) is void*, and ArrayFire does not export array::device<void*>(),
// so device<CUdeviceptr>() would be an unresolved external. Fetch the raw device pointer through an
// exported element type (float) and reinterpret it to CUdeviceptr* (void**); the pointer value is the
// device address, matching how the CUDA path uses device<CUdeviceptr>() (CUdeviceptr is an integer there).
#define transferAF(varA) reinterpret_cast<CUdeviceptr*>((varA).device<float>())
#else
#define transferAF(varA) (varA).device<CUdeviceptr>()
#endif
#elif defined(METAL)
#include "ProjectorClass.h"
#define transferAF(varA) NS::RetainPtr(reinterpret_cast<MTL::Buffer*>((varA).device<void>()))
#elif defined(CPU)
#include "ProjectorClassCPU.h"
#define transferAF(varA) varA.device<float>()
#else
#include "ProjectorClass.h"
#define transferAF(varA) cl::Buffer(*varA.device<cl_mem>(), true)
#endif
#define VAL 0.00001f
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#pragma warning(disable : 4996)

/// <summary>
/// Determines whether the image-domain part of the subset-based PDHG/PKMA/MBSREM/BSREM update can be fused into the
/// projector type 4/5 backprojection kernel. The fused kernel never allocates the backprojection and skips the
/// algorithm specific functions and regularization after backprojection (except in the case of BSREM)
/// The fused step is used automatically whenever the configuration supports it; the decision is silent unless
/// DEBUG is defined, in which case a message states whether it was enabled and, if not, why. An explicit
/// inputScalars.fastPDHG = false disables it regardless of whether the configuration would otherwise support it.
/// </summary>
/// <param name="inputScalars various scalar parameters. inputScalars.fastPDHG is updated in place"></param>
/// <param name="w_vec specifies algorithm and prior parameters"></param>
/// <param name="MethodList specifies the algorithms and priors used"></param>
inline void checkFastPDHG(scalarStruct& inputScalars, const Weighting& w_vec, const RecMethods& MethodList) {
	if (!inputScalars.fastPDHG) {
		// Explicitly disabled by the user
		if (DEBUG)
			mexPrint("fastPDHG disabled: Turned off explicitly by the user.");
		return;
	}
#ifdef CPU
	inputScalars.fastPDHG = false;
#elif defined(METAL)
	// The fused step has no Metal implementation
	inputScalars.fastPDHG = false;
	if (DEBUG)
		mexPrint("fastPDHG disabled: not supported by the Metal backend.");
	return;
#else
	const char* reason = nullptr;
	// fastPDHG supports a limited number of algorithms: PDHG (and its KL/L1/CV variants) and some
	// Poisson-based algorithms PKMA/MBSREM/BSREM
	const bool fastPDHGAlgorithm = (MethodList.PDHG || MethodList.PDHGKL || MethodList.PDHGL1 || MethodList.CV) && !MethodList.PDDY;
	const bool fastPoissonAlgorithm = MethodList.PKMA || MethodList.MBSREM || MethodList.BSREM;
	// Only the voxel-driven type 4 and type 5 backprojections compute one voxel per work-item, which is what
	// makes the image-domain update possible in the first place
	// All BPType 4/5 combinations are allowed
	if ((inputScalars.BPType != 4 && inputScalars.BPType != 5) || !inputScalars.CT)
		reason = "only backprojection types 4 and 5 (projector_type 4, 5, 14, 15, 24, 34, 45, 54, 74, 75) with CT data are supported";
	else if (!fastPDHGAlgorithm && !fastPoissonAlgorithm)
		reason = "only PDHG, PDHGKL, PDHGL1, CV, PKMA, MBSREM and BSREM are supported (PDDY is not)";
	// The experimental relaxation parameter computations are not supported
	else if (MethodList.PKMA && (inputScalars.computeRelaxation || inputScalars.relaxScaling))
		reason = "PKMA with computeRelaxation or relaxScaling is not supported by the fused step";
	// fastPDHG is aimed at CT and thus there's no support for listmode
	else if (fastPoissonAlgorithm && inputScalars.listmode > 0)
		reason = "listmode data is not supported with the fused PKMA/MBSREM/BSREM step";
	else if (inputScalars.subsets <= 1 || inputScalars.stochastic)
		reason = "only the subset-based case is supported";
	else if (inputScalars.Nt > 1)
		reason = "dynamic (multiple time step) data is not supported";
	else if (inputScalars.adaptiveType > 0)
		reason = "adaptive step sizes (PDAdaptiveType) are not supported";
	else if (inputScalars.use_psf)
		reason = "PSF is not supported";
	// No atomics should be used with the voxel-based BPs
	else if (inputScalars.atomic_64bit || inputScalars.atomic_32bit)
		reason = "integer atomics are not supported";
	else if (inputScalars.FISTAAcceleration)
		reason = "FISTA acceleration is not supported";
	// When using extended FOV, the regularization is computed only for the main FOV/volume region
	// This is currently not supported with fastPDHG, however, multi-resolution can be used as the
	// volumes are computed separately
	else if ((inputScalars.useExtendedFOV || inputScalars.eFOV) && !inputScalars.multiResolution)
		reason = "the extended FOV without multi-resolution is not supported";
	// Only a limited number of regularization methods are supported, mainly due to the need to have a
	// separate kernel
	// Proximal TV/TGV are not included due to the more complicated structure despite having separate kernels
	else if (MethodList.MRP || MethodList.Quad || MethodList.Huber || MethodList.L || MethodList.FMH || MethodList.AD ||
		MethodList.APLS || MethodList.TGV || MethodList.WeightedMean ||
		MethodList.ProxTV || MethodList.ProxTGV || MethodList.ProxRDP || MethodList.ProxNLM)
		reason = "the only supported priors are NLM, RDP, GGMRF, hyperbolic and TV";
	// Anatomical reference images are not supported in any case
	// TODO: Add anatomical reference support? Not that relevant in CT
	else if (MethodList.TV && (w_vec.data.TV_use_anatomical || w_vec.data.TVtype == 6))
		reason = "anatomical/weighted TV is not supported by the fused step";
	// The use of images/textures is important with fastPDHG so those need to be enabled
	else if ((MethodList.NLM || MethodList.RDP || MethodList.GGMRF || MethodList.hyperbolic || MethodList.TV) && !inputScalars.useImages)
		reason = "NLM/RDP/GGMRF/hyperbolic/TV with the fused PDHG step requires images (useImages = true)";
	else if (MethodList.TemporalSmoothness || MethodList.TemporalTV)
		reason = "temporal priors are not supported";
	// Preconditioners 0 and 1 are computed in the kernel, the rest are not
	else if (w_vec.precondTypeIm[2] || w_vec.precondTypeIm[3] || w_vec.precondTypeIm[4] || w_vec.precondTypeIm[5] || w_vec.precondTypeIm[6])
		reason = "only the image-based preconditioners 0 and 1 are supported";
	else if (w_vec.precondTypeIm[0] && w_vec.precondTypeIm[1])
		reason = "only one image-based preconditioner can be used at a time";
	if (reason != nullptr) {
		inputScalars.fastPDHG = false;
		// The fused step is enabled automatically whenever it applies, so a configuration that does not
		// support it is not an error, print only for debugging purposes
		if (DEBUG) {
			mexPrintBase("fastPDHG disabled: %s. Using the regular (multi-step) path.\n", reason);
			mexEval();
		}
	}
	else {
		if (DEBUG)
			mexPrint("fastPDHG enabled: the image-domain update is fused into the backprojection kernel");
	}
#endif
}

/// <summary>
/// Sets the per-sub-iteration parameters needed by fastPDHG backprojection kernel
/// </summary>
#ifndef CPU
static inline void setFastPDHGSubIterParams(ProjectorClass& proj, const scalarStruct& inputScalars, const Weighting& w_vec,
	const RecMethods& MethodList, const int tt, const int iter, const int osa_iter) {
	const int64_t curIter = static_cast<int64_t>(osa_iter) + static_cast<int64_t>(inputScalars.subsets) * static_cast<int64_t>(iter);
	const int64_t curIterMax = static_cast<int64_t>(inputScalars.subsets) * static_cast<int64_t>(inputScalars.Niter) - 1;
	const bool computeReg = inputScalars.regEveryIter <= 1 || curIter == 0 || curIter == curIterMax
		|| (curIter % static_cast<int64_t>(inputScalars.regEveryIter)) == 0;
	proj.fastStep = 1;
	// PDHG
	if (proj.fastAlg == 0) {
		proj.fastTheta = static_cast<float>(inputScalars.subsetsUsed) * w_vec.thetaCP[tt][curIter];
		proj.fastTau = w_vec.tauCP[tt][0];
	}
	else {
		// Poisson algorithms (PKMA/MBSREM/BSREM)
		// Relaxation parameter
		proj.fastLambda = w_vec.lambda[tt][iter];
		if (MethodList.PKMA) {
			// Momentum parameter
			const uint32_t kk = iter * inputScalars.subsets + inputScalars.currentSubset;
			proj.fastAlpha = w_vec.alphaM[tt][kk];
		}
		else if (MethodList.MBSREM)
			// Upper bound
			proj.fastAlpha = w_vec.U;
		else // BSREM
			proj.fastAlpha = 1.f;
	}
	proj.fastEpps = inputScalars.epps;
	proj.fastPositivity = inputScalars.enforcePositivity ? 1 : 0;
	// If the regularization is not computed, set beta as zero
	// BSREM doesn't compute regularization as part of the backprojection kernel
	proj.fastBeta = (proj.fastNLMUsed && computeReg && !MethodList.BSREM) ? w_vec.beta : 0.f;
}
#endif

// This function sets the variables needed for the special large dimensional reconstruction method
inline void largeDimCreate(scalarStruct& inputScalars, const RecMethods& MethodList, const Weighting& w_vec) {
	inputScalars.lDimStruct.Nz.resize(inputScalars.subsets);
	inputScalars.lDimStruct.NzPr.resize(inputScalars.subsets);
	inputScalars.lDimStruct.startPr.resize(inputScalars.subsets);
	inputScalars.lDimStruct.endPr.resize(inputScalars.subsets);
	inputScalars.lDimStruct.imDim.resize(inputScalars.subsets);
	inputScalars.lDimStruct.imDimPr.resize(inputScalars.subsets);
	inputScalars.lDimStruct.cumDim.resize(inputScalars.subsets + 1);
	inputScalars.lDimStruct.cumDimPr.resize(inputScalars.subsets + 1);
	inputScalars.lDimStruct.bz.resize(inputScalars.subsets);
	inputScalars.lDimStruct.bmaxZ.resize(inputScalars.subsets);
	inputScalars.lDimStruct.d_Scale4Z.resize(inputScalars.subsets);
	inputScalars.lDimStruct.cumDim[0] = 0LL;
	inputScalars.lDimStruct.cumDimPr[0] = 0LL;
	const uint32_t intZ = inputScalars.Nz[0] / inputScalars.subsets;
	const uint32_t remainder = inputScalars.Nz[0] % inputScalars.subsets;
	for (int kk = 0; kk < inputScalars.subsets; kk++) {
		if (kk == 0) {
			inputScalars.lDimStruct.Nz[kk] = intZ + remainder;
			inputScalars.lDimStruct.bz[kk] = inputScalars.bz[0];
		}
		else {
			inputScalars.lDimStruct.Nz[kk] = intZ;
			inputScalars.lDimStruct.bz[kk] = inputScalars.lDimStruct.bmaxZ[kk - 1];
		}
		inputScalars.lDimStruct.bmaxZ[kk] = inputScalars.lDimStruct.bz[kk] + inputScalars.lDimStruct.Nz[kk] * inputScalars.dz[0];
		if (MethodList.NLM || MethodList.RDP || MethodList.hyperbolic || MethodList.GGMRF || MethodList.TV) {
			inputScalars.lDimStruct.imDim[kk] = static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]) * static_cast<int64_t>(inputScalars.lDimStruct.Nz[kk]);
			inputScalars.lDimStruct.cumDim[kk + 1] = inputScalars.lDimStruct.cumDim[kk] + inputScalars.lDimStruct.imDim[kk];
			if (kk < inputScalars.subsetsUsed - 1 && kk > 0) {
				if (MethodList.NLM) {
					inputScalars.lDimStruct.imDimPr[kk] = static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]) * (static_cast<int64_t>(inputScalars.lDimStruct.Nz[kk] + w_vec.Ndz * 2U + w_vec.Nlz * 2U));
					inputScalars.lDimStruct.NzPr[kk] = inputScalars.lDimStruct.Nz[kk] + w_vec.Ndz * 2U + w_vec.Nlz * 2U;
				}
				else if (MethodList.TV || (MethodList.RDP && !w_vec.RDPLargeNeighbor)) {
					inputScalars.lDimStruct.imDimPr[kk] = static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]) * (static_cast<int64_t>(inputScalars.lDimStruct.Nz[kk] + 2U));
					inputScalars.lDimStruct.NzPr[kk] = inputScalars.lDimStruct.Nz[kk] + 2U;
				}
				else {
					inputScalars.lDimStruct.imDimPr[kk] = static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]) * (static_cast<int64_t>(inputScalars.lDimStruct.Nz[kk] + w_vec.Ndz * 2U));
					inputScalars.lDimStruct.NzPr[kk] = inputScalars.lDimStruct.Nz[kk] + w_vec.Ndz * 2U;
				}
			}
			else {
				if (MethodList.NLM) {
					inputScalars.lDimStruct.imDimPr[kk] = static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]) * (static_cast<int64_t>(inputScalars.lDimStruct.Nz[kk] + w_vec.Ndz + w_vec.Nlz));
					inputScalars.lDimStruct.NzPr[kk] = inputScalars.lDimStruct.Nz[kk] + w_vec.Ndz + w_vec.Nlz;
				}
				else if (MethodList.TV || (MethodList.RDP && !w_vec.RDPLargeNeighbor)) {
					inputScalars.lDimStruct.imDimPr[kk] = static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]) * (static_cast<int64_t>(inputScalars.lDimStruct.Nz[kk] + 1));
					inputScalars.lDimStruct.NzPr[kk] = inputScalars.lDimStruct.Nz[kk] + 1;
				}
				else {
					inputScalars.lDimStruct.imDimPr[kk] = static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]) * (static_cast<int64_t>(inputScalars.lDimStruct.Nz[kk] + w_vec.Ndz));
					inputScalars.lDimStruct.NzPr[kk] = inputScalars.lDimStruct.Nz[kk] + w_vec.Ndz;
				}
			}
			if (kk == 0) {
				if (MethodList.NLM)
					inputScalars.lDimStruct.cumDimPr[kk + 1] = inputScalars.lDimStruct.imDim[kk] - static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]) * (w_vec.Ndz + w_vec.Nlz);
				else if (MethodList.TV || (MethodList.RDP && !w_vec.RDPLargeNeighbor))
					inputScalars.lDimStruct.cumDimPr[kk + 1] = inputScalars.lDimStruct.imDim[kk] - static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]);
				else
					inputScalars.lDimStruct.cumDimPr[kk + 1] = inputScalars.lDimStruct.imDim[kk] - static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]) * w_vec.Ndz;
				inputScalars.lDimStruct.startPr[kk] = 0;
			}
			else
				if (MethodList.NLM) {
					inputScalars.lDimStruct.cumDimPr[kk + 1] = inputScalars.lDimStruct.cumDim[kk] + inputScalars.lDimStruct.imDim[kk] - static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]) * (w_vec.Ndz + w_vec.Nlz);
					inputScalars.lDimStruct.startPr[kk] = (w_vec.Ndz + w_vec.Nlz) * static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]);
				}
				else if (MethodList.TV || (MethodList.RDP && !w_vec.RDPLargeNeighbor)) {
					inputScalars.lDimStruct.startPr[kk] = static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]);
					inputScalars.lDimStruct.cumDimPr[kk + 1] = inputScalars.lDimStruct.cumDim[kk] + inputScalars.lDimStruct.imDim[kk] - static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]);
				}
				else {
					inputScalars.lDimStruct.startPr[kk] = w_vec.Ndz * static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]);
					inputScalars.lDimStruct.cumDimPr[kk + 1] = inputScalars.lDimStruct.cumDim[kk] + inputScalars.lDimStruct.imDim[kk] - static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]) * w_vec.Ndz;
				}
			if (kk == inputScalars.subsetsUsed - 1)
				inputScalars.lDimStruct.endPr[kk] = 0U;
			else
				if (MethodList.NLM)
					inputScalars.lDimStruct.endPr[kk] = (w_vec.Ndz + w_vec.Nlz) * static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]);
				else if (MethodList.TV || (MethodList.RDP && !w_vec.RDPLargeNeighbor))
					inputScalars.lDimStruct.endPr[kk] = static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]);
				else
					inputScalars.lDimStruct.endPr[kk] = w_vec.Ndz * static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]);
		}
		else {
			inputScalars.lDimStruct.imDim[kk] = static_cast<int64_t>(inputScalars.Nx[0]) * static_cast<int64_t>(inputScalars.Ny[0]) * static_cast<int64_t>(inputScalars.lDimStruct.Nz[kk]);
			inputScalars.lDimStruct.cumDim[kk + 1] = inputScalars.lDimStruct.cumDim[kk] + inputScalars.lDimStruct.imDim[kk];
			inputScalars.lDimStruct.imDimPr[kk] = inputScalars.lDimStruct.imDim[kk];
			inputScalars.lDimStruct.cumDimPr[kk + 1] = inputScalars.lDimStruct.cumDim[kk + 1];
			inputScalars.lDimStruct.NzPr[kk] = inputScalars.lDimStruct.Nz[kk];
		}
		inputScalars.lDimStruct.d_Scale4Z[kk] = 1.f / (static_cast<float>(inputScalars.lDimStruct.Nz[kk]) * inputScalars.dz[0]);
		if (DEBUG) {
			mexPrintBase("inputScalars.lDimStruct.NzPr[kk] = %u\n", inputScalars.lDimStruct.NzPr[kk]);
			mexPrintBase("inputScalars.lDimStruct.Nz[kk] = %u\n", inputScalars.lDimStruct.Nz[kk]);
			mexPrintBase("inputScalars.lDimStruct.imDim[kk] = %u\n", inputScalars.lDimStruct.imDim[kk]);
			mexPrintBase("inputScalars.lDimStruct.imDimPr[kk] = %u\n", inputScalars.lDimStruct.imDimPr[kk]);
			mexEval();
		}
	}
}

inline void largeDimFirst(scalarStruct& inputScalars, ProjectorClass& proj, const int iter) {
	if (iter == 0) {
		inputScalars.lDimStruct.NzOrig = inputScalars.Nz[0];
		inputScalars.lDimStruct.imDimOrig = inputScalars.im_dim[0];
#ifndef CPU
		inputScalars.lDimStruct.bzOrig = VEC_Z(proj.b[0]);
		inputScalars.lDimStruct.bmaxZOrig = VEC_Z(proj.bmax[0]);
		inputScalars.lDimStruct.d_Scale4ZOrig = VEC_Z(inputScalars.d_Scale4[0]);
#endif
	}
	inputScalars.Nz[0] = inputScalars.lDimStruct.Nz[iter];
	inputScalars.im_dim[0] = inputScalars.lDimStruct.imDim[iter];
#ifndef CPU
	VEC_Z(proj.d_N[0]) = inputScalars.lDimStruct.Nz[iter];
	VEC_Z(proj.b[0]) = inputScalars.lDimStruct.bz[iter];
	VEC_Z(proj.bmax[0]) = inputScalars.lDimStruct.bmaxZ[iter];
	VEC_Z(inputScalars.d_Scale4[0]) = inputScalars.lDimStruct.d_Scale4Z[iter];
#endif
}

inline void largeDimLast(scalarStruct& inputScalars, ProjectorClass& proj) {
	inputScalars.Nz[0] = inputScalars.lDimStruct.NzOrig;
	inputScalars.im_dim[0] = inputScalars.lDimStruct.imDimOrig;
#ifndef CPU
	VEC_Z(proj.d_N[0]) = inputScalars.lDimStruct.NzOrig;
	VEC_Z(proj.b[0]) = inputScalars.lDimStruct.bzOrig;
	VEC_Z(proj.bmax[0]) = inputScalars.lDimStruct.bmaxZOrig;
	VEC_Z(inputScalars.d_Scale4[0]) = inputScalars.lDimStruct.d_Scale4ZOrig;
#endif
}

// Fill the backprojection (aka right-hand side) with zeros
inline void initializeRHS(AF_im_vectors& vec, const scalarStruct& inputScalars, uint32_t timestep, const int ii = 0) {
	if (inputScalars.verbose >= 3)
		mexPrint("Initialize the backprojection output");
#if defined(OPENCL) || defined(METAL)
	if (inputScalars.atomic_64bit)
		vec.rhs_os[timestep][ii] = (af::constant(0LL, static_cast<size_t>(inputScalars.im_dim[ii]) * static_cast<size_t>(inputScalars.nRekos), 1, s64));
	else if (inputScalars.atomic_32bit)
		vec.rhs_os[timestep][ii] = (af::constant(0, static_cast<size_t>(inputScalars.im_dim[ii]) * static_cast<size_t>(inputScalars.nRekos), 1, s32));
	else
#endif
		vec.rhs_os[timestep][ii] = (af::constant(0.f, static_cast<size_t>(inputScalars.im_dim[ii]) * static_cast<size_t>(inputScalars.nRekos), 1));
	vec.rhs_os[timestep][ii].eval();
}

// The backprojection mask of the SPECT projector type 6 with multi-resolution volumes, for every volume (w_vec.maskBPVol[frame][volume], 0 or 1, the
// voxels of the volume with x the fastest), made once on the host. The mask is one mask of the whole (total) grid of the reconstruction, which applies
// to every volume: the union bounding box of the volumes at the voxel size of volume 0 (the grid of options.x0 and of the attenuation map, see
// buildType6AttenuationVolume), 2D (maskBPZ == 1, the same for every slice) or 3D. A voxel of a volume is inside the mask when its centre, in world
// coordinates, is in a voxel of the whole-grid mask that is inside: a centre on the boundary of two voxels of the grid takes the one with the higher
// index (as the nearest-neighbour interpolation of MATLAB does) and a centre outside the grid is outside the mask. The multi-resolution host code of
// the MATLAB interface of OMEGA can give the masks already as one mask for each volume (the volumes one after the other, always 3D, and one such set
// for every frame when there is a mask for every frame of a dynamic scan; maskBPZ is then the largest number of slices of the volumes), which is used
// as it is: that layout is recognised by the number of elements. When the number of elements is not known (size_maskBP 0) the mask is taken to be
// of the whole grid. Does nothing without multi-resolution volumes (the mask is the one of the volume, see maskBPArray). Returns false (after
// printing an error) if the mask has neither layout
inline bool prepareMaskBPType6(const scalarStruct& inputScalars, Weighting& w_vec) {
	w_vec.maskBPVol.clear();
	if (w_vec.maskBP == nullptr || inputScalars.nMultiVolumes == 0)
		return true;
	const int nV = static_cast<int>(inputScalars.nMultiVolumes) + 1;
	std::vector<size_t> vox(static_cast<size_t>(nV));
	size_t sumVox = 0ULL;
	uint32_t maxNz = 0U;
	for (int ii = 0; ii < nV; ii++) {
		vox[ii] = static_cast<size_t>(inputScalars.Nx[ii]) * static_cast<size_t>(inputScalars.Ny[ii]) * static_cast<size_t>(inputScalars.Nz[ii]);
		sumVox += vox[ii];
		maxNz = (std::max)(maxNz, static_cast<uint32_t>(inputScalars.Nz[ii]));
	}
	const size_t sizeMask = inputScalars.size_maskBP;
	if (sizeMask >= sumVox && sizeMask % sumVox == 0ULL && inputScalars.maskBPZ > 1U && inputScalars.maskBPZ == maxNz &&
		(sizeMask / sumVox == 1ULL || sizeMask / sumVox == static_cast<size_t>(inputScalars.Nt))) {
		// One mask for each volume
		const size_t frames = sizeMask / sumVox;
		w_vec.maskBPVol.assign(frames, std::vector<std::vector<uint8_t>>(static_cast<size_t>(nV)));
		size_t off = 0ULL;
		for (size_t ff = 0; ff < frames; ff++) {
			for (int ii = 0; ii < nV; ii++) {
				w_vec.maskBPVol[ff][ii].assign(w_vec.maskBP + off, w_vec.maskBP + off + vox[ii]);
				off += vox[ii];
			}
		}
		return true;
	}
	// One mask of the whole grid
	float bxMap = inputScalars.bx[0], byMap = inputScalars.by[0], bzMap = inputScalars.bz[0];
	float bxHi = inputScalars.bx[0] + inputScalars.Nx[0] * inputScalars.dx[0];
	float byHi = inputScalars.by[0] + inputScalars.Ny[0] * inputScalars.dy[0];
	float bzHi = inputScalars.bz[0] + inputScalars.Nz[0] * inputScalars.dz[0];
	for (int jj = 1; jj < nV; jj++) {
		bxMap = (std::min)(bxMap, inputScalars.bx[jj]);
		byMap = (std::min)(byMap, inputScalars.by[jj]);
		bzMap = (std::min)(bzMap, inputScalars.bz[jj]);
		bxHi = (std::max)(bxHi, inputScalars.bx[jj] + inputScalars.Nx[jj] * inputScalars.dx[jj]);
		byHi = (std::max)(byHi, inputScalars.by[jj] + inputScalars.Ny[jj] * inputScalars.dy[jj]);
		bzHi = (std::max)(bzHi, inputScalars.bz[jj] + inputScalars.Nz[jj] * inputScalars.dz[jj]);
	}
	const int NxMap = static_cast<int>(std::lround((bxHi - bxMap) / inputScalars.dx[0]));
	const int NyMap = static_cast<int>(std::lround((byHi - byMap) / inputScalars.dy[0]));
	const int NzMap = static_cast<int>(std::lround((bzHi - bzMap) / inputScalars.dz[0]));
	const bool mask3D = inputScalars.maskBPZ > 1U;
	const size_t need = static_cast<size_t>(NxMap) * static_cast<size_t>(NyMap) * (mask3D ? static_cast<size_t>(NzMap) : 1ULL);
	if (sizeMask != 0ULL && sizeMask != need) {
		mexPrintBase("maskBP: with multi-resolution volumes the mask must be of the whole grid, %d x %d", NxMap, NyMap);
		if (mask3D)
			mexPrintBase(" x %d", NzMap);
		mexPrintBase(" voxels, but it has %llu elements\n", static_cast<unsigned long long>(sizeMask));
		mexEval();
		return false;
	}
	// The voxel of the whole grid that the centre of every voxel of the volume is in (-1: outside the grid), per axis. The small offset puts the centres
	// that are on the boundary of two voxels of the grid (a double size voxel on the grid of the voxels of the main volume) in the voxel with the higher index
	const auto gridIndex = [](const double b, const double d, const double bMap, const double dMap, const int n, const int nMap) {
		std::vector<int> idx(static_cast<size_t>(n));
		for (int i = 0; i < n; i++) {
			const double u = (b + (static_cast<double>(i) + 0.5) * d - bMap) / dMap;
			const long long k = static_cast<long long>(std::floor(u + 1e-4));
			idx[i] = (k >= 0LL && k < static_cast<long long>(nMap)) ? static_cast<int>(k) : -1;
		}
		return idx;
	};
	std::vector<std::vector<uint8_t>> vols(static_cast<size_t>(nV));
	for (int ii = 0; ii < nV; ii++) {
		const int nx = static_cast<int>(inputScalars.Nx[ii]), ny = static_cast<int>(inputScalars.Ny[ii]), nz = static_cast<int>(inputScalars.Nz[ii]);
		const std::vector<int> ix = gridIndex(inputScalars.bx[ii], inputScalars.dx[ii], bxMap, inputScalars.dx[0], nx, NxMap);
		const std::vector<int> iy = gridIndex(inputScalars.by[ii], inputScalars.dy[ii], byMap, inputScalars.dy[0], ny, NyMap);
		const std::vector<int> iz = gridIndex(inputScalars.bz[ii], inputScalars.dz[ii], bzMap, inputScalars.dz[0], nz, NzMap);
		vols[ii].assign(vox[ii], 0);
		size_t inside = 0ULL;
		for (int kk = 0; kk < nz; kk++) {
			if (mask3D && iz[kk] < 0)
				continue;
			const size_t zOff = mask3D ? static_cast<size_t>(iz[kk]) * static_cast<size_t>(NxMap) * static_cast<size_t>(NyMap) : 0ULL;
			for (int jj = 0; jj < ny; jj++) {
				if (iy[jj] < 0)
					continue;
				for (int i = 0; i < nx; i++) {
					if (ix[i] < 0)
						continue;
					if (w_vec.maskBP[zOff + static_cast<size_t>(iy[jj]) * static_cast<size_t>(NxMap) + static_cast<size_t>(ix[i])] != 0) {
						vols[ii][static_cast<size_t>(i) + static_cast<size_t>(nx) * (static_cast<size_t>(jj) + static_cast<size_t>(ny) * static_cast<size_t>(kk))] = 1;
						inside++;
					}
				}
			}
		}
		if (DEBUG || inputScalars.verbose >= 3)
			mexPrintBase("Backprojection mask of volume %d: %llu of %llu voxels inside\n", ii, static_cast<unsigned long long>(inside), static_cast<unsigned long long>(vox[ii]));
	}
	w_vec.maskBPVol.assign(1, vols);
	return true;
}

// Multiplies a flattened (Nx*Ny*Nz) volume by the backprojection mask. The mask is either 2D (Nx x Ny, maskBPZ == 1, applied to every slice)
// or fully 3D (Nx x Ny x maskBPZ). Returns false (after printing an error) if the host mask is too small, to avoid reading out of bounds.
// With the masks of the volumes of prepareMaskBPType6 (the type-6 projector with multi-resolution volumes) the mask of the volume ii (of the
// frame timestep, when there is a mask for every frame) is that one
inline bool maskBPArray(af::array& maskBP, const scalarStruct& inputScalars, const Weighting& w_vec, const int ii, const int timestep = 0) {
	if (!w_vec.maskBPVol.empty()) {
		const size_t nFrames = w_vec.maskBPVol.size();
		const size_t frame = nFrames > 1U ? (std::min)(static_cast<size_t>((std::max)(timestep, 0)), nFrames - 1U) : static_cast<size_t>(0);
		const std::vector<uint8_t>& vol = w_vec.maskBPVol[frame][static_cast<size_t>(ii)];
		maskBP = af::array(inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii], vol.data());
		maskBP = af::flat((maskBP != 0).as(f32));
		return true;
	}
	const size_t Nx = static_cast<size_t>(inputScalars.Nx[ii]), Ny = static_cast<size_t>(inputScalars.Ny[ii]), Nz = static_cast<size_t>(inputScalars.Nz[ii]);
	const bool mask3D = inputScalars.maskBPZ > 1U;
	const size_t need = mask3D ? Nx * Ny * Nz : Nx * Ny;
	if (inputScalars.size_maskBP != 0ULL && inputScalars.size_maskBP < need) {
		mexPrintBase("maskBP: host array has %llu elements but ", static_cast<unsigned long long>(inputScalars.size_maskBP));
		mexPrintBase("%llu are required\n", static_cast<unsigned long long>(need));
		mexEval();
		return false;
	}
	if (mask3D)
		maskBP = af::array(inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii], w_vec.maskBP);
	else
		maskBP = af::tile(af::array(inputScalars.Nx[ii], inputScalars.Ny[ii], w_vec.maskBP), 1, 1, inputScalars.Nz[ii]);
	maskBP = af::flat((maskBP != 0).as(f32));
	return true;
}

inline bool applyMaskBP(af::array& output, af::array* outputSens, const scalarStruct& inputScalars, const Weighting& w_vec, const int ii, const int timestep = 0) {
	af::array maskBP;
	if (!maskBPArray(maskBP, inputScalars, w_vec, ii, timestep))
		return false;
	output *= maskBP;
	if (outputSens != nullptr)
		*outputSens *= maskBP;
	return true;
}

// Gives the voxels outside the backprojection mask the values that the projectors of types 1-3 leave there after their backprojection and sensitivity
// finalization: nothing in the measurement backprojection (rhs = 0, never raised to epps) and the smallest sensitivity (epps). The multiplicative update
// im / Summ * (rhs + epps) of OSEM then leaves such voxels at exactly their current value. Call this after the epps clamps of the type-6 backprojection
// (which would otherwise turn the zero of the masked rhs into epps and the zero of the masked sensitivity into 1)
inline bool resetMaskedBP(af::array& rhs, af::array* Summ, const scalarStruct& inputScalars, const Weighting& w_vec, const int ii, const int timestep = 0) {
	af::array maskBP;
	if (!maskBPArray(maskBP, inputScalars, w_vec, ii, timestep))
		return false;
	const af::array inside = maskBP != 0.f;
	rhs = af::select(inside, rhs, 0.f);
	if (Summ != nullptr)
		*Summ = af::select(inside, *Summ, inputScalars.epps);
	return true;
}

// Zero and symmetric padding of an array
inline af::array padding(const af::array& im, const uint32_t Nx, const uint32_t Ny, const uint32_t Nz, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz,
	const bool zero_pad = false, const uint32_t Nw = 1)
{
	af::array out;
	if (zero_pad == 1) {
		af::dtype type = im.type();
		if (Nz == 1) {
			//if (in.dims(1) == 1)
				//in = af::moddims(in, Nx, Ny, Nz, Nw);
			out = af::constant(0, Nx + 2 * Ndx, Ny + 2 * Ndy, 1, Nw, type);
			out(static_cast<double>(Ndx) + af::seq(static_cast<double>(out.dims(0))), static_cast<double>(Ndy) + af::seq(static_cast<double>(out.dims(1))), af::span, af::span) = af::moddims(im, Nx, Ny, Nz, Nw);
		}
		else {
			//if (out.dims(2) == 1)
				//out = af::moddims(out, Nx, Ny, Nz, Nw);
			out = af::constant(0, Nx + 2 * Ndx, Ny + 2 * Ndy, Nz + 2 * Ndz, Nw, type);
			out(static_cast<double>(Ndx) + af::seq(static_cast<double>(out.dims(0))), static_cast<double>(Ndy) + af::seq(static_cast<double>(out.dims(1))),
				static_cast<double>(Ndz) + af::seq(static_cast<double>(out.dims(2))), af::span) = af::moddims(im, Nx, Ny, Nz, Nw);
		}
	}
	else {
		out = im;
		if (out.dims(1) == 1)
			out = af::moddims(out, Nx, Ny, Nz, Nw);
		// Symmetric (mirrored) padding
		// The manual flip/join version built the padding with af::flip, but that fails on POCL devices
		// Requires AF 3.8 or newer
		const dim_t padX = static_cast<dim_t>(Ndx);
		const dim_t padY = static_cast<dim_t>(Ndy);
		const dim_t padZ = (Nz == 1) ? 0LL : static_cast<dim_t>(Ndz);
		if (padX > 0 || padY > 0 || padZ > 0)
			out = af::pad(out, af::dim4(padX, padY, padZ, 0), af::dim4(padX, padY, padZ, 0), AF_PAD_SYM);
	}
	return out;
}

// Compute the convolution of an array
inline af::array computeConvolution(const af::array& vec, const af::array& g, const scalarStruct& inputScalars, const Weighting& w_vec,
	const uint32_t nRekos = 1, const int ii = 0) {
	if (inputScalars.verbose >= 3)
		mexPrint("Starting PSF blurring");
	af::array apu = af::moddims(vec, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii], nRekos);
	apu = padding(apu, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii], inputScalars.g_dim_x + 1, inputScalars.g_dim_y + 1, inputScalars.g_dim_z + 1,
		false, inputScalars.nRekos2);
	apu = af::convolve3(apu, g);
	if (inputScalars.Nz[ii] > 1)
		apu = apu(af::seq(inputScalars.g_dim_x + 1, inputScalars.Nx[ii] + inputScalars.g_dim_x), af::seq(inputScalars.g_dim_y + 1, inputScalars.Ny[ii] + inputScalars.g_dim_y),
			af::seq(inputScalars.g_dim_z + 1, inputScalars.Nz[ii] + inputScalars.g_dim_z), af::span);
	else
		apu = apu(af::seq(inputScalars.g_dim_x + 1, inputScalars.Nx[ii] + inputScalars.g_dim_x), af::seq(inputScalars.g_dim_y + 1, inputScalars.Ny[ii] + inputScalars.g_dim_y),
			af::span, af::span);
	if (inputScalars.verbose >= 3)
		mexPrint("PSF blurring complete");
	return af::flat(apu);
}

/// <summary>
/// Transfer the sensitivity image from an ArrayFire array to an OpenCL buffer or CUDA device pointer
/// </summary>
/// <param name="apuSum the sensitivity image"></param>
/// <param name="proj the projector class object"></param>
/// <returns></returns>
inline void transferSensitivityImage(af::array& apuSum, ProjectorClass& proj) {
	apuSum.eval();
	if (proj.d_Summ.size() < 1)
		proj.d_Summ.emplace_back(transferAF(apuSum));
	else
		proj.d_Summ[0] = transferAF(apuSum);
}

/// <summary>
/// Transfer the backprojection from an ArrayFire array to an OpenCL buffer or CUDA device pointer
/// </summary>
/// <param name="rhs_os the backprojection"></param>
/// <param name="proj the projector class object"></param>
/// <returns></returns>
inline int transferRHS(af::array& rhs_os, ProjectorClass& proj) {
	if (DEBUG) {
		mexPrintBase("proj.vec_opencl.d_rhs_os.size() = %u\n", proj.vec_opencl.d_rhs_os.size());
		mexEval();
	}
	if (proj.vec_opencl.d_rhs_os.size() < 1)
		proj.vec_opencl.d_rhs_os.emplace_back(transferAF(rhs_os));
	else
		proj.vec_opencl.d_rhs_os[0] = transferAF(rhs_os);
	if (DEBUG) {
		mexPrintBase("proj.vec_opencl.d_rhs_os.size() = %u\n", proj.vec_opencl.d_rhs_os.size());
		mexEval();
	}
	return 0;
}

/// <summary>
/// Copy the current estimates from an ArrayFire array to an OpenCL image/CUDA texture. For branchless distance-driven, compute the integral images before copy
/// </summary>
/// <param name="vec image estimates and backprojection"></param>
/// <param name="inputScalars various scalar parameters defining the build parameters and what features to use"></param>
/// <param name="proj the projector class object"></param>
/// <param name="ii optional multi-resolution volume number, default is 0 (main volume)"></param>
/// <returns></returns>
inline int updateInputs(AF_im_vectors& vec, const scalarStruct& inputScalars, ProjectorClass& proj, const int ii, const int timestep) {
#if defined(CUDA) || defined(HIP)
	CUresult status = CUDA_SUCCESS;
#elif defined(OPENCL)
	int status = 0;
	cl::detail::size_t_array region = { inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii] };
#elif defined(METAL)
	int status = 0;
#endif
	if (inputScalars.FPType == 5) {
		af::array im;
		af::array intIm = af::constant(0.f, inputScalars.Ny[ii] + 1, inputScalars.Nz[ii] + 1, inputScalars.Nx[ii]);
		if (inputScalars.meanFP) {
			im = af::reorder(af::moddims(vec.im_os[timestep][ii], inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii], inputScalars.nRekos), 1, 2, 0);
			const af::array mFPXY = af::mean(af::mean(im, 0), 1);
			vec.meanFP = af::constant(0.f, inputScalars.Nx[ii] + inputScalars.Ny[ii]);
			vec.meanFP(af::seq(0, inputScalars.Nx[ii] - 1)) = af::flat(mFPXY);
			im -= af::tile(mFPXY, im.dims(0), im.dims(1));
			intIm(af::seq(1, af::end), af::seq(1, af::end), af::span) = af::sat(im);
			im.eval();
		}
		else
			intIm(af::seq(1, af::end), af::seq(1, af::end), af::span) = af::sat(af::reorder(af::moddims(vec.im_os[timestep][ii], inputScalars.Nx[ii],
				inputScalars.Ny[ii], inputScalars.Nz[ii] * inputScalars.nRekos), 1, 2, 0));
		intIm.eval();
		dim_t dim0 = intIm.dims(0);
		dim_t dim1 = intIm.dims(1);
		dim_t dim2 = intIm.dims(2);
		if (DEBUG) {
			mexPrintBase("dim0 = %u\n", dim0);
			mexPrintBase("dim1 = %u\n", dim1);
			mexPrintBase("dim2 = %u\n", dim2);
			mexPrintBase("af::sum<float>(intIm) = %f\n", af::sum<float>(intIm));
			mexPrintBase("af::sum<float>(vec.im_os[timestep][ii]) = %f\n", af::sum<float>(vec.im_os[timestep][ii]));
			mexEval();
		}
		intIm = af::flat(intIm);
#if defined(CUDA) || defined(HIP)
		if (proj.intTexCacheXY.size() <= static_cast<size_t>(ii)) {
			proj.intTexCacheXY.resize(ii + 1);
			proj.intArrayCacheXY.resize(ii + 1);
			proj.intImageCacheDimsXY.resize((ii + 1) * 3, 0);
		}
		if (proj.intImageCacheDimsXY[ii * 3] != static_cast<size_t>(dim0) || proj.intImageCacheDimsXY[ii * 3 + 1] != static_cast<size_t>(dim1) || proj.intImageCacheDimsXY[ii * 3 + 2] != static_cast<size_t>(dim2)) {
			if (proj.intImageCacheDimsXY[ii * 3 + 2] != 0) {
				getErrorString(cuTexObjectDestroy(proj.intTexCacheXY[ii]));
				getErrorString(cuArrayDestroy(proj.intArrayCacheXY[ii]));
			}
			CUDA_TEXTURE_DESC texDesc;
			CUDA_ARRAY3D_DESCRIPTOR_st arr3DDesc;
			CUDA_RESOURCE_DESC resDesc;
			CUDA_RESOURCE_VIEW_DESC viewDesc;
			std::memset(&texDesc, 0, sizeof(texDesc));
			std::memset(&resDesc, 0, sizeof(resDesc));
			std::memset(&arr3DDesc, 0, sizeof(arr3DDesc));
			std::memset(&viewDesc, 0, sizeof(viewDesc));
			arr3DDesc.Format = CUarray_format::CU_AD_FORMAT_FLOAT;
			arr3DDesc.NumChannels = 1;
			arr3DDesc.Height = dim1;
			arr3DDesc.Width = dim0;
			arr3DDesc.Depth = dim2;
			status = cuArray3DCreate(&proj.intArrayCacheXY[ii], &arr3DDesc);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create integral array\n");
				return -1;
			}
			else if (DEBUG)
				mexPrint("Arrray creation completed\n");
			resDesc.resType = CUresourcetype::CU_RESOURCE_TYPE_ARRAY;
			resDesc.res.array.hArray = proj.intArrayCacheXY[ii];
			texDesc.addressMode[0] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
			texDesc.addressMode[1] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
			texDesc.addressMode[2] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
			texDesc.filterMode = CUfilter_mode::CU_TR_FILTER_MODE_LINEAR;
			texDesc.flags = CU_TRSF_NORMALIZED_COORDINATES;
			viewDesc.height = dim1;
			viewDesc.width = dim0;
			viewDesc.depth = dim2;
			viewDesc.format = CUresourceViewFormat::CU_RES_VIEW_FORMAT_FLOAT_1X32;
			status = cuTexObjectCreate(&proj.intTexCacheXY[ii], &resDesc, &texDesc, &viewDesc);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				mexPrint("Integral image xz copy failed\n");
				return -1;
			}
			else if (DEBUG)
				mexPrint("Texture creation completed\n");
			proj.intImageCacheDimsXY[ii * 3] = static_cast<size_t>(dim0);
			proj.intImageCacheDimsXY[ii * 3 + 1] = static_cast<size_t>(dim1);
			proj.intImageCacheDimsXY[ii * 3 + 2] = static_cast<size_t>(dim2);
		}
		proj.vec_opencl.d_image_os_int = proj.intTexCacheXY[ii];
		CUDA_MEMCPY3D cpy3d;
		std::memset(&cpy3d, 0, sizeof(cpy3d));
		cpy3d.srcMemoryType = CUmemorytype::CU_MEMORYTYPE_DEVICE;
		cpy3d.srcDevice = reinterpret_cast<CUdeviceptr>(transferAF(intIm));
		cpy3d.srcPitch = dim0 * sizeof(float);
		cpy3d.srcHeight = dim1;
		cpy3d.dstMemoryType = CUmemorytype::CU_MEMORYTYPE_ARRAY;
		cpy3d.dstArray = proj.intArrayCacheXY[ii];
		cpy3d.WidthInBytes = dim0 * sizeof(float);
		cpy3d.Height = dim1;
		cpy3d.Depth = dim2;
		status = cuMemcpy3DAsync(&cpy3d, proj.CLCommandQueue[0]);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to copy integral array\n");
			return -1;
		}
		else if (DEBUG)
			mexPrint("Arrray copy completed\n");
#elif defined(OPENCL)
		region = { static_cast<cl_ulong>(dim0), static_cast<cl_ulong>(dim1), static_cast<cl_ulong>(dim2) };
		if (proj.d_intImageCacheXY.size() <= static_cast<size_t>(ii)) {
			proj.d_intImageCacheXY.resize(ii + 1);
			proj.intImageCacheDimsXY.resize((ii + 1) * 3, 0);
		}
		if (proj.intImageCacheDimsXY[ii * 3] != region[0] || proj.intImageCacheDimsXY[ii * 3 + 1] != region[1] || proj.intImageCacheDimsXY[ii * 3 + 2] != region[2]) {
			proj.d_intImageCacheXY[ii] = cl::Image3D(proj.CLContext, CL_MEM_READ_ONLY, proj.format, region[0], region[1], region[2], 0, 0, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create integral image\n");
				return -1;
			}
			proj.intImageCacheDimsXY[ii * 3] = region[0];
			proj.intImageCacheDimsXY[ii * 3 + 1] = region[1];
			proj.intImageCacheDimsXY[ii * 3 + 2] = region[2];
		}
		proj.vec_opencl.d_image_os_int = proj.d_intImageCacheXY[ii];
		cl_int status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(cl::Buffer(*intIm.device<cl_mem>(), true), proj.vec_opencl.d_image_os_int, 0, proj.origin, region);
		if (status != 0) {
			getErrorString(status);
			mexPrint("Integral image xz copy failed\n");
			return -1;
		}
#elif defined(METAL)
		{
			auto source = transferAF(intIm);
			status = proj.makeImageTextureFromDevice(proj.vec_opencl.d_image_os_int, proj.imArray, source,
				static_cast<size_t>(dim0), static_cast<size_t>(dim1), static_cast<size_t>(dim2),
				BACKEND_TEXTURE_LINEAR, BACKEND_TEXTURE_NORMALIZED);
			if (status != 0) {
				mexPrint("Integral image xz copy failed\n");
				return -1;
			}
		}
#endif
		intIm.unlock();
		intIm = af::constant(0.f, inputScalars.Nx[ii] + 1, inputScalars.Nz[ii] + 1, inputScalars.Ny[ii]);
		if (inputScalars.meanFP) {
			im = af::reorder(af::moddims(vec.im_os[timestep][ii], inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii], inputScalars.nRekos), 0, 2, 1, 3);
			const af::array mFPYZ = af::mean(af::mean(im, 0), 1);
			vec.meanFP(af::seq(inputScalars.Nx[ii], inputScalars.Nx[ii] + inputScalars.Ny[ii] - 1)) = af::flat(mFPYZ);
			im -= af::tile(mFPYZ, im.dims(0), im.dims(1));
			intIm(af::seq(1, af::end), af::seq(1, af::end), af::span) = af::sat(im);
		}
		else
			intIm(af::seq(1, af::end), af::seq(1, af::end), af::span) = af::sat(af::reorder(af::moddims(vec.im_os[timestep][ii], inputScalars.Nx[ii],
				inputScalars.Ny[ii], inputScalars.Nz[ii], inputScalars.nRekos), 0, 2, 1, 3));
		dim0 = intIm.dims(0);
		dim1 = intIm.dims(1);
		dim2 = intIm.dims(2);
		intIm = af::flat(intIm);
#if defined(CUDA) || defined(HIP)
		if (proj.intTexCacheYZ.size() <= static_cast<size_t>(ii)) {
			proj.intTexCacheYZ.resize(ii + 1);
			proj.intArrayCacheYZ.resize(ii + 1);
			proj.intImageCacheDimsYZ.resize((ii + 1) * 3, 0);
		}
		if (proj.intImageCacheDimsYZ[ii * 3] != static_cast<size_t>(dim0) || proj.intImageCacheDimsYZ[ii * 3 + 1] != static_cast<size_t>(dim1) || proj.intImageCacheDimsYZ[ii * 3 + 2] != static_cast<size_t>(dim2)) {
			if (proj.intImageCacheDimsYZ[ii * 3 + 2] != 0) {
				getErrorString(cuTexObjectDestroy(proj.intTexCacheYZ[ii]));
				getErrorString(cuArrayDestroy(proj.intArrayCacheYZ[ii]));
			}
			CUDA_TEXTURE_DESC texDesc;
			CUDA_ARRAY3D_DESCRIPTOR_st arr3DDesc;
			CUDA_RESOURCE_DESC resDesc;
			CUDA_RESOURCE_VIEW_DESC viewDesc;
			std::memset(&texDesc, 0, sizeof(texDesc));
			std::memset(&resDesc, 0, sizeof(resDesc));
			std::memset(&arr3DDesc, 0, sizeof(arr3DDesc));
			std::memset(&viewDesc, 0, sizeof(viewDesc));
			arr3DDesc.Format = CUarray_format::CU_AD_FORMAT_FLOAT;
			arr3DDesc.NumChannels = 1;
			arr3DDesc.Height = dim1;
			arr3DDesc.Width = dim0;
			arr3DDesc.Depth = dim2;
			status = cuArray3DCreate(&proj.intArrayCacheYZ[ii], &arr3DDesc);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create integral array\n");
				return -1;
			}
			else if (DEBUG)
				mexPrint("Arrray creation completed\n");
			resDesc.resType = CUresourcetype::CU_RESOURCE_TYPE_ARRAY;
			resDesc.res.array.hArray = proj.intArrayCacheYZ[ii];
			texDesc.addressMode[0] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
			texDesc.addressMode[1] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
			texDesc.addressMode[2] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
			texDesc.filterMode = CUfilter_mode::CU_TR_FILTER_MODE_LINEAR;
			texDesc.flags = CU_TRSF_NORMALIZED_COORDINATES;
			viewDesc.height = dim1;
			viewDesc.width = dim0;
			viewDesc.depth = dim2;
			viewDesc.format = CUresourceViewFormat::CU_RES_VIEW_FORMAT_FLOAT_1X32;
			status = cuTexObjectCreate(&proj.intTexCacheYZ[ii], &resDesc, &texDesc, &viewDesc);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				mexPrint("Integral image yz copy failed\n");
				return -1;
			}
			else if (DEBUG)
				mexPrint("Texture creation completed\n");
			proj.intImageCacheDimsYZ[ii * 3] = static_cast<size_t>(dim0);
			proj.intImageCacheDimsYZ[ii * 3 + 1] = static_cast<size_t>(dim1);
			proj.intImageCacheDimsYZ[ii * 3 + 2] = static_cast<size_t>(dim2);
		}
		proj.vec_opencl.d_image_os = proj.intTexCacheYZ[ii];
		CUDA_MEMCPY3D cpy3dYZ;
		std::memset(&cpy3dYZ, 0, sizeof(cpy3dYZ));
		cpy3dYZ.srcMemoryType = CUmemorytype::CU_MEMORYTYPE_DEVICE;
		cpy3dYZ.srcDevice = reinterpret_cast<CUdeviceptr>(transferAF(intIm));
		cpy3dYZ.srcPitch = dim0 * sizeof(float);
		cpy3dYZ.srcHeight = dim1;
		cpy3dYZ.dstMemoryType = CUmemorytype::CU_MEMORYTYPE_ARRAY;
		cpy3dYZ.dstArray = proj.intArrayCacheYZ[ii];
		cpy3dYZ.WidthInBytes = dim0 * sizeof(float);
		cpy3dYZ.Height = dim1;
		cpy3dYZ.Depth = dim2;
		status = cuMemcpy3DAsync(&cpy3dYZ, proj.CLCommandQueue[0]);
		if (status != CUDA_SUCCESS) {
			getErrorString(status);
			mexPrint("Failed to copy integral array\n");
			return -1;
		}
		else if (DEBUG)
			mexPrint("Arrray copy completed\n");
#elif defined(OPENCL)
		region = { static_cast<cl_ulong>(dim0), static_cast<cl_ulong>(dim1), static_cast<cl_ulong>(dim2) };
		if (proj.d_intImageCacheYZ.size() <= static_cast<size_t>(ii)) {
			proj.d_intImageCacheYZ.resize(ii + 1);
			proj.intImageCacheDimsYZ.resize((ii + 1) * 3, 0);
		}
		if (proj.intImageCacheDimsYZ[ii * 3] != region[0] || proj.intImageCacheDimsYZ[ii * 3 + 1] != region[1] || proj.intImageCacheDimsYZ[ii * 3 + 2] != region[2]) {
			proj.d_intImageCacheYZ[ii] = cl::Image3D(proj.CLContext, CL_MEM_READ_ONLY, proj.format, region[0], region[1], region[2], 0, 0, NULL, &status);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to create integral image\n");
				return -1;
			}
			proj.intImageCacheDimsYZ[ii * 3] = region[0];
			proj.intImageCacheDimsYZ[ii * 3 + 1] = region[1];
			proj.intImageCacheDimsYZ[ii * 3 + 2] = region[2];
		}
		proj.vec_opencl.d_image_os = proj.d_intImageCacheYZ[ii];
		status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(cl::Buffer(*intIm.device<cl_mem>(), true), proj.vec_opencl.d_image_os, 0, proj.origin, region);
		if (status != 0) {
			getErrorString(status);
			mexPrint("Integral image yz copy failed\n");
			return -1;
		}
#elif defined(METAL)
		{
			auto source = transferAF(intIm);
			status = proj.makeImageTextureFromDevice(proj.vec_opencl.d_image_os, proj.imArray, source,
				static_cast<size_t>(dim0), static_cast<size_t>(dim1), static_cast<size_t>(dim2),
				BACKEND_TEXTURE_LINEAR, BACKEND_TEXTURE_NORMALIZED);
			if (status != 0) {
				mexPrint("Integral image yz copy failed\n");
				return -1;
			}
		}
#endif
		intIm.unlock();
	}
	else {
#if defined(CUDA) || defined(HIP)
		if (inputScalars.useBuffers) {
			if (inputScalars.use_psf)
				proj.vec_opencl.d_im = transferAF(vec.im_os_blurred[ii]);
			else
				proj.vec_opencl.d_im = transferAF(vec.im_os[timestep][ii]);
		}
		else {
			CUdeviceptr* im;
			if (inputScalars.use_psf)
				im = transferAF(vec.im_os_blurred[ii]);
			else
				im = transferAF(vec.im_os[timestep][ii]);
			if (proj.FPTexCache.size() <= static_cast<size_t>(ii)) {
				proj.FPTexCache.resize(ii + 1);
				proj.FPArrayCache.resize(ii + 1);
				proj.FPTexCachePrior.resize(ii + 1);
				proj.imageCacheDims.resize((ii + 1) * 3, 0);
			}
			if (proj.imageCacheDims[ii * 3] != static_cast<size_t>(inputScalars.Nx[ii]) || proj.imageCacheDims[ii * 3 + 1] != static_cast<size_t>(inputScalars.Ny[ii]) || proj.imageCacheDims[ii * 3 + 2] != static_cast<size_t>(inputScalars.Nz[ii])) {
				if (proj.imageCacheDims[ii * 3 + 2] != 0) {
					getErrorString(cuTexObjectDestroy(proj.FPTexCache[ii]));
					getErrorString(cuTexObjectDestroy(proj.FPTexCachePrior[ii]));
					getErrorString(cuArrayDestroy(proj.FPArrayCache[ii]));
				}
				CUDA_TEXTURE_DESC texDesc;
				CUDA_ARRAY3D_DESCRIPTOR_st arr3DDesc;
				CUDA_RESOURCE_DESC resDesc;
				CUDA_RESOURCE_VIEW_DESC viewDesc;
				std::memset(&texDesc, 0, sizeof(texDesc));
				std::memset(&resDesc, 0, sizeof(resDesc));
				std::memset(&arr3DDesc, 0, sizeof(arr3DDesc));
				std::memset(&viewDesc, 0, sizeof(viewDesc));
				arr3DDesc.Format = CUarray_format::CU_AD_FORMAT_FLOAT;
				arr3DDesc.NumChannels = 1;
				arr3DDesc.Height = inputScalars.Ny[ii];
				arr3DDesc.Width = inputScalars.Nx[ii];
				arr3DDesc.Depth = inputScalars.Nz[ii];
				status = cuArray3DCreate(&proj.FPArrayCache[ii], &arr3DDesc);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					mexPrint("Failed to create image array\n");
					return -1;
				}
				resDesc.resType = CUresourcetype::CU_RESOURCE_TYPE_ARRAY;
				resDesc.res.array.hArray = proj.FPArrayCache[ii];
				if (inputScalars.FPType == 4) {
					texDesc.addressMode[0] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
					texDesc.addressMode[1] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
					texDesc.addressMode[2] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
					texDesc.filterMode = CUfilter_mode::CU_TR_FILTER_MODE_LINEAR;
					texDesc.flags = CU_TRSF_NORMALIZED_COORDINATES;
				}
				else {
					texDesc.addressMode[0] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
					texDesc.addressMode[1] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
					texDesc.addressMode[2] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
					texDesc.filterMode = CUfilter_mode::CU_TR_FILTER_MODE_POINT;
				}
				viewDesc.height = inputScalars.Ny[ii];
				viewDesc.width = inputScalars.Nx[ii];
				viewDesc.depth = inputScalars.Nz[ii];
				viewDesc.format = CUresourceViewFormat::CU_RES_VIEW_FORMAT_FLOAT_1X32;
				status = cuTexObjectCreate(&proj.FPTexCache[ii], &resDesc, &texDesc, &viewDesc);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					mexPrint("Image texture creation failed\n");
					status = cuArrayDestroy(proj.FPArrayCache[ii]);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
					}
					return -1;
				}
				// If FPType == 4, the texture fetching uses normalized coordinates while the priors require integers ones
				// Create another description of the same texture using integer coordinates with nearest neighbor interpolation
				CUDA_TEXTURE_DESC texDescPrior;
				std::memset(&texDescPrior, 0, sizeof(texDescPrior));
				texDescPrior.addressMode[0] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
				texDescPrior.addressMode[1] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
				texDescPrior.addressMode[2] = CUaddress_mode::CU_TR_ADDRESS_MODE_CLAMP;
				texDescPrior.filterMode = CUfilter_mode::CU_TR_FILTER_MODE_POINT;
				status = cuTexObjectCreate(&proj.FPTexCachePrior[ii], &resDesc, &texDescPrior, &viewDesc);
				if (status != CUDA_SUCCESS) {
					getErrorString(status);
					mexPrint("Prior image texture creation failed\n");
					getErrorString(cuTexObjectDestroy(proj.FPTexCache[ii]));
					status = cuArrayDestroy(proj.FPArrayCache[ii]);
					if (status != CUDA_SUCCESS) {
						getErrorString(status);
					}
					return -1;
				}
				proj.imageCacheDims[ii * 3] = static_cast<size_t>(inputScalars.Nx[ii]);
				proj.imageCacheDims[ii * 3 + 1] = static_cast<size_t>(inputScalars.Ny[ii]);
				proj.imageCacheDims[ii * 3 + 2] = static_cast<size_t>(inputScalars.Nz[ii]);
			}
			proj.vec_opencl.d_image_os = proj.FPTexCache[ii];
			if (DEBUG) {
				mexPrintBase("vec.im_os[timestep][ii].elements() = %u\n", vec.im_os[timestep][ii].elements());
				mexPrintBase("inputScalars.Nx[ii] * sizeof(float) = %u\n", inputScalars.Nx[ii] * sizeof(float));
				mexPrintBase("inputScalars.Ny[ii] = %u\n", inputScalars.Ny[ii]);
				mexPrintBase("inputScalars.Nz[ii] = %u\n", inputScalars.Nz[ii]);
				mexPrintBase("elements = %u\n", inputScalars.Nz[ii] * inputScalars.Nx[ii] * inputScalars.Ny[ii]);
				mexEval();
			}
			CUDA_MEMCPY3D cpy3d;
			std::memset(&cpy3d, 0, sizeof(cpy3d));
			cpy3d.srcMemoryType = CUmemorytype::CU_MEMORYTYPE_DEVICE;
			cpy3d.srcDevice = reinterpret_cast<CUdeviceptr>(im);
			cpy3d.srcPitch = inputScalars.Nx[ii] * sizeof(float);
			cpy3d.srcHeight = inputScalars.Ny[ii];
			cpy3d.dstMemoryType = CUmemorytype::CU_MEMORYTYPE_ARRAY;
			cpy3d.dstArray = proj.FPArrayCache[ii];
			cpy3d.WidthInBytes = inputScalars.Nx[ii] * sizeof(float);
			cpy3d.Height = inputScalars.Ny[ii];
			cpy3d.Depth = inputScalars.Nz[ii];
			status = cuMemcpy3DAsync(&cpy3d, proj.CLCommandQueue[0]);
			if (status != CUDA_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to copy image array\n");
				return -1;
			}
			if (inputScalars.use_psf)
				vec.im_os_blurred[ii].unlock();
			else
				vec.im_os[timestep][ii].unlock();
		}
#elif defined(OPENCL)
		if (inputScalars.useBuffers) {
			if (inputScalars.use_psf)
				proj.vec_opencl.d_im = cl::Buffer(*vec.im_os_blurred[ii].device<cl_mem>(), true);
			else
				proj.vec_opencl.d_im = cl::Buffer(*vec.im_os[timestep][ii].device<cl_mem>(), true);
		}
		else {
			if (proj.d_imageCache.size() <= static_cast<size_t>(ii)) {
				proj.d_imageCache.resize(ii + 1);
				proj.imageCacheDims.resize((ii + 1) * 3, 0);
			}
			if (proj.imageCacheDims[ii * 3] != region[0] || proj.imageCacheDims[ii * 3 + 1] != region[1] || proj.imageCacheDims[ii * 3 + 2] != region[2]) {
				proj.d_imageCache[ii] = cl::Image3D(proj.CLContext, CL_MEM_READ_ONLY, proj.format, region[0], region[1], region[2], 0, 0, NULL, &status);
				if (status != CL_SUCCESS) {
					getErrorString(status);
					mexPrint("Failed to create input images\n");
					return -1;
				}
				else if (DEBUG)
					mexPrint("Input image created\n");
				proj.imageCacheDims[ii * 3] = region[0];
				proj.imageCacheDims[ii * 3 + 1] = region[1];
				proj.imageCacheDims[ii * 3 + 2] = region[2];
			}
			proj.vec_opencl.d_image_os = proj.d_imageCache[ii];
			cl_mem* im;
			if (inputScalars.use_psf)
				im = vec.im_os_blurred[ii].device<cl_mem>();
			else
				im = vec.im_os[timestep][ii].device<cl_mem>();
			cl::Buffer d_im_os = cl::Buffer(*im, true);
			status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(d_im_os, proj.vec_opencl.d_image_os, 0, proj.origin, region);
			// Mod: no finish needed, same in-order queue as ArrayFire
			if (inputScalars.use_psf)
				vec.im_os_blurred[ii].unlock();
			else
				vec.im_os[timestep][ii].unlock();
			delete im;
			if (status != 0) {
				getErrorString(status);
				mexPrint("Image copy failed\n");
				return -1;
			}
			else if (DEBUG)
				mexPrint("Input copy succeeded\n");
		}
#elif defined(METAL)
		if (inputScalars.use_psf)
			proj.vec_opencl.d_im = transferAF(vec.im_os_blurred[ii]);
		else
			proj.vec_opencl.d_im = transferAF(vec.im_os[timestep][ii]);
		proj.vec_opencl.d_image_os.reset(); // The texture wraps the current volume buffer and dimensions.
#else
		if (inputScalars.use_psf)
			proj.vec_opencl.d_im_os = vec.im_os_blurred[ii].device<float>();
		else
			proj.vec_opencl.d_im_os = vec.im_os[timestep][ii].device<float>();
#endif
	}
	return 0;
}

/// <summary>
/// Computes the square root of the offset redundancy weights at the detector pixel centers for every subset. The weights match
/// the ones computed in the projector type 4/5 backprojection kernels when OFFSETSQRT is defined, i.e. sqrt(.5 * (1 + sin(pi * px / (2 * T))))
/// </summary>
/// <param name="inputScalars various scalar parameters, inputScalars.T contains the distance (mm) of the central ray from the row-0 edge of the detector for each projection"></param>
/// <param name="w_vec the weights are stored to w_vec.offsetSqrtW (nRowsD x 1 x number of projections for each subset)"></param>
/// <param name="length number of projections in each subset"></param>
/// <param name="pituus cumulative number of projections"></param>
inline void computeOffsetSqrtWeights(const scalarStruct& inputScalars, Weighting& w_vec, const std::vector<int64_t>& length, const int64_t* pituus) {
	w_vec.offsetSqrtW.assign(inputScalars.subsetsUsed, af::array());
	const float koko = static_cast<float>(inputScalars.nRowsD) * w_vec.dPitchX;
	for (uint32_t kk = inputScalars.osa_iter0; kk < inputScalars.subsetsUsed; kk++) {
		const dim_t nP = length[kk];
		const af::array T = af::tile(af::array(1, 1, nP, &inputScalars.T[pituus[kk]], afHost), inputScalars.nRowsD, 1, 1);
		// Normalized detector coordinate of the pixel centers
		af::array px = (af::range(af::dim4(inputScalars.nRowsD, 1, nP), 0, f32) + 0.5f) / static_cast<float>(inputScalars.nRowsD);
		// Same flipping as in the kernels when the central ray is on the far half of the detector
		const af::array flip = T > koko * 0.5f;
		px = af::select(flip, af::abs(px - 1.f), px);
		const af::array TT = af::select(flip, koko - T, T);
		px = px * koko - TT;
		af::array w = af::sin(3.14159265358979f * (px + TT) / (4.f * TT));
		w = af::select(px > TT, 1.f, w);
		w = af::select(px < -TT, 0.f, w);
		w.eval();
		w_vec.offsetSqrtW[kk] = w;
	}
}

/// <summary>
/// Multiplies a measurement-domain array (nRowsD x nColsD x projections of the subset kk) by the square root of the offset redundancy weights
/// </summary>
/// <param name="m the array, replaced by the weighted (single precision) array"></param>
/// <param name="inputScalars various scalar parameters"></param>
/// <param name="w_vec contains the weights computed with computeOffsetSqrtWeights"></param>
/// <param name="kk the subset"></param>
inline void applyOffsetSqrtW(af::array& m, const scalarStruct& inputScalars, const Weighting& w_vec, const uint32_t kk) {
	const af::array& W = w_vec.offsetSqrtW[kk];
	const af::dim4 dd = m.dims();
	m = af::moddims(af::moddims(m.as(f32), inputScalars.nRowsD, inputScalars.nColsD, W.dims(2)) * af::tile(W, 1, inputScalars.nColsD, 1), dd);
	m.eval();
}

// Helper function to transfer AF arrays into device vectors, call the actual forward projection function and move memory management back to AF
// If inputScalars.offsetSqrt is enabled, the forward projection is weighted with the square root of the offset redundancy weights (the backprojection kernels apply
// the same weights on the other side). rawFP, if given, receives the unweighted forward projection (added to it)
inline int forwardProjectionAFOpenCL(AF_im_vectors& vec, scalarStruct& inputScalars, Weighting& w_vec, af::array& outputFP, uint32_t osa_iter, uint32_t timestep,
	const std::vector<int64_t>& length, const af::array& g, uint64_t m_size, ProjectorClass& proj, const int ii = 0, const int64_t* pituus = nullptr, af::array* rawFP = nullptr) {
	int status = 0;
	if (inputScalars.use_psf) {
		const int nRekos = vec.im_os[timestep][ii].elements() / (inputScalars.im_dim[ii]);
		vec.im_os_blurred[ii] = computeConvolution(vec.im_os[timestep][ii], g, inputScalars, w_vec, nRekos, ii);
	}
	if (DEBUG) {
		mexPrintBase("outputFP.dims(0) = %d\n", outputFP.dims(0));
		mexPrintBase("outputFP.dims(1) = %d\n", outputFP.dims(1));
		mexEval();
	}
	// The forward projection kernels accumulate into the output (and with multiple volumes every volume adds to the same output), so the weighting is applied
	// to a zero-initialized temporary that is then added to the output
	const bool sqrtFP = inputScalars.offsetSqrt && osa_iter < w_vec.offsetSqrtW.size() && !w_vec.offsetSqrtW[osa_iter].isempty();
	af::array fpTemp;
	if (sqrtFP) {
		if (static_cast<dim_t>(outputFP.elements()) != static_cast<dim_t>(inputScalars.nRowsD) * static_cast<dim_t>(inputScalars.nColsD) * w_vec.offsetSqrtW[osa_iter].dims(2)) {
			mexPrint("The forward projection size does not match the offset weights");
			return -1;
		}
		fpTemp = af::constant(0.f, outputFP.dims(), f32);
	}
	af::array& fpTarget = sqrtFP ? fpTemp : outputFP;
	proj.d_output = transferAF(fpTarget);
	status = updateInputs(vec, inputScalars, proj, ii, timestep);
	if (status != 0) {
		return -1;
	}
	proj.memSize += (sizeof(float) * inputScalars.im_dim[ii]);
#ifndef CPU
	if (inputScalars.meanFP && inputScalars.FPType == 5)
		proj.d_meanFP = transferAF(vec.meanFP);
	status = proj.forwardProjection(inputScalars, w_vec, osa_iter, timestep, length, m_size, ii);
	if (inputScalars.useBuffers) {
		if (inputScalars.use_psf)
			vec.im_os_blurred[ii].unlock();
		else
			vec.im_os[timestep][ii].unlock();
	}
#else
	status = proj.forwardProjection(inputScalars, w_vec, osa_iter, timestep, length, pituus, ii);
	if (inputScalars.use_psf)
		vec.im_os_blurred[ii].unlock();
	else
		vec.im_os[timestep][ii].unlock();
#endif
	fpTarget.unlock();
	if (inputScalars.meanFP && inputScalars.FPType == 5)
		vec.meanFP.unlock();
	proj.memSize -= (sizeof(float) * inputScalars.im_dim[ii]);
	if (sqrtFP && status == 0) {
		if (rawFP != nullptr)
			*rawFP += fpTemp;
		applyOffsetSqrtW(fpTemp, inputScalars, w_vec, osa_iter);
		outputFP += fpTemp;
		outputFP.eval();
	}
	return status;
}

// Same as above, but for backprojection
inline int backwardProjectionAFOpenCL(AF_im_vectors& vec, scalarStruct& inputScalars, Weighting& w_vec, const RecMethods& MethodList, af::array& outputFP, uint32_t osa_iter, uint32_t timestep,
	std::vector<int64_t>& length, uint64_t m_size, af::array& meanBP, const af::array& g, ProjectorClass& proj, const bool compSens = false, const int ii = 0,
	const int64_t* pituus = nullptr, const bool FDK = false, const int queueIdx = 0, const bool finalize = true, const bool newInput = true, const int lDimChunk = -1) {
	int status = 0;
	outputFP.eval();
	// The fastPDHG directly updates the image estimate, so it needs the estimate itself, (possibly) the PDHG auxiliary image
	// and, when applicable, the sensitivity image of the preconditioner and the NLM patch weights
    // When fastPDHG is used, it writes the new estimate directly into im_os, so the usual RHS (backprojection output)
	// initialization is skipped
#ifndef CPU
	const bool fastStep = inputScalars.fastPDHG && proj.fastStep != 0;
	const bool fastPrecond = w_vec.precondTypeIm[0] || w_vec.precondTypeIm[1];
#else
	// fastPDHG is not supported by the CPU implementation
	const bool fastStep = false;
#endif
	if (!FDK && !fastStep)
		initializeRHS(vec, inputScalars, timestep, ii);
	proj.memSize += (sizeof(float) * inputScalars.im_dim[ii]);
	if (DEBUG) {
		mexPrintBase("ii = %u\n", ii);
		mexPrintBase("vec.rhs_os[timestep][ii].dims(0) = %u\n", vec.rhs_os[timestep][ii].dims(0));
		mexPrintBase("inputScalars.nRekos2 = %u\n", inputScalars.nRekos2);
		mexPrintBase("inputScalars.nRekos = %u\n", inputScalars.nRekos);
		mexPrintBase("outputFP = %f\n", af::sum<float>(outputFP));
		mexPrintBase("min(outputFP) = %f\n", af::min<float>(outputFP));
		mexPrintBase("max(outputFP) = %f\n", af::max<float>(outputFP));
		mexEval();
	}
	proj.d_output = transferAF(outputFP);
	if (DEBUG) {
		mexPrint("Transferring backprojection output\n");
	}
	if (fastStep)
		status = transferRHS(vec.im_os[timestep][ii], proj);
	else
		status = transferRHS(vec.rhs_os[timestep][ii], proj);
	if (status != 0) {
		return -1;
	}
	if (DEBUG)
		mexPrint("Backprojection output transfered\n");
#ifndef CPU
	// Only transfer if NLM is used. The patch weights are a single shared array (not per volume) and only the
	// main volume is regularized, so lock them once, for ii == 0. Locking them again for every volume would
	// leave several locks to be released by the single unlock below (concurrent multi-volume case).
	if (ii == 0 && inputScalars.fastPDHG && proj.fastNLMUsed && MethodList.NLM)
		proj.d_gaussianNLM = transferAF(w_vec.gaussianNLM);
	if (fastStep) {
		// There is no uCP (PDHG backprojection accumulator) for the Poisson algorithms (PKMA/MBSREM/BSREM)
		if (proj.fastAlg == 0)
			proj.d_U = transferAF(vec.uCP[timestep][ii]);
		// Primal step-size is per volume
		if (proj.fastAlg == 0)
			proj.fastTau = w_vec.tauCP[timestep][ii];
		// Preconditioner vector (sensitivity image)
		if (fastPrecond)
			proj.d_precond = transferAF(w_vec.D[0][ii]);
		// Include the subvolume information of largeDim
		if (inputScalars.largeDim && lDimChunk >= 0) {
			proj.fastImOffset = static_cast<int64_t>(inputScalars.lDimStruct.startPr[lDimChunk]);
			proj.fastPriorZOffset = static_cast<int>(inputScalars.lDimStruct.startPr[lDimChunk] / (inputScalars.Nx[0] * inputScalars.Ny[0]));
		}
		else {
			proj.fastImOffset = 0;
			proj.fastPriorZOffset = 0;
		}
		// Priors supported by fastPDHG require a copy of the current estimate as an image/texture
		// When using FPType == 5, the same image cannot be used as BDD uses integral images
		// In such cases, the image needs to be copied beforehand
		// This also applies to largeDim cases where the estimate is computed with subvolumes/chunks
		// Only the main volume is regularized, so the estimate image is created once, for ii == 0
		if (ii == 0 && proj.fastNLMUsed && (inputScalars.FPType == 5 || inputScalars.largeDim)) {
#if defined(CUDA) || defined(HIP)
			uint32_t Nz = inputScalars.Nz[0];
			if (inputScalars.largeDim && lDimChunk >= 0)
				Nz = inputScalars.lDimStruct.NzPr[lDimChunk];
			if (proj.transferTex(inputScalars, proj.vec_opencl.d_rhs_os[0], false, Nz) != 0)
				return -1;
#elif defined(OPENCL)
			const size_t origRegionZ = proj.region[2];
			if (inputScalars.largeDim && lDimChunk >= 0) {
				proj.region[2] = inputScalars.lDimStruct.NzPr[lDimChunk];
				proj.d_inputI = cl::Image3D(proj.CLContext, CL_MEM_READ_ONLY, proj.format, proj.region[0], proj.region[1], proj.region[2], 0, 0, NULL, &status);
				OCL_CHECK(status, "Failed to create prior image\n", -1);
			}
			status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(proj.vec_opencl.d_rhs_os[0], proj.d_inputI, 0, proj.origin, proj.region);
			if (status != CL_SUCCESS) {
				getErrorString(status);
				mexPrint("Failed to copy the estimate into the prior image\n");
				return -1;
			}
			if (inputScalars.largeDim && lDimChunk >= 0)
				proj.region[2] = origRegionZ;
#endif
		}
	}
	if (inputScalars.meanBP && inputScalars.BPType == 5)
		proj.d_meanBP = transferAF(meanBP);
#if defined(CUDA) || defined(HIP)
	status = proj.backwardProjection(inputScalars, w_vec, osa_iter, timestep, length, m_size, MethodList, compSens, ii, 0, queueIdx, newInput);
#elif defined(OPENCL) || defined(METAL)
	status = proj.backwardProjection(inputScalars, w_vec, osa_iter, timestep, length, m_size, MethodList, compSens, ii, 0, -1, queueIdx, newInput);
#endif
	// Only unlock after all the queues/streams have finalized
	if (fastStep && finalize) {
		if (proj.fastAlg == 0)
			vec.uCP[timestep][ii].unlock();
		vec.im_os[timestep][ii].unlock();
		if (fastPrecond)
			w_vec.D[0][ii].unlock();
#if defined(CUDA) || defined(HIP)
		// Destroy the texture if it's created above (only ever created for ii == 0)
		// TODO: Persistent image/texture for fastPDHG?
		if (ii == 0 && proj.fastNLMUsed && (inputScalars.FPType == 5 || inputScalars.largeDim)) {
			CUresult statusTex = cuTexObjectDestroy(proj.d_inputI);
			if (statusTex != CUDA_SUCCESS)
				getErrorString(statusTex);
			statusTex = cuArrayDestroy(proj.imArray);
			if (statusTex != CUDA_SUCCESS)
				getErrorString(statusTex);
		}
#endif
	}
	if (ii == 0 && finalize && inputScalars.fastPDHG && proj.fastNLMUsed && MethodList.NLM)
		w_vec.gaussianNLM.unlock();
#else
	 status = proj.backwardProjection(inputScalars, w_vec, osa_iter, timestep, length, pituus, compSens, ii);
#endif
	if (!finalize)
		return status;
	if (!fastStep)
		vec.rhs_os[timestep][ii].unlock();
	outputFP.unlock();
	if (inputScalars.meanBP && inputScalars.BPType == 5)
		meanBP.unlock();
	if (status != 0) {
		return -1;
	}
	if (!fastStep) {
#if defined(OPENCL) || defined(METAL)
		if (inputScalars.atomic_64bit)
			vec.rhs_os[timestep][ii] = vec.rhs_os[timestep][ii].as(f32) / TH;
		else if (inputScalars.atomic_32bit)
			vec.rhs_os[timestep][ii] = vec.rhs_os[timestep][ii].as(f32) / TH32;
#endif
		if (inputScalars.use_psf) {
			const int nRekos = vec.rhs_os[timestep][ii].elements() / (inputScalars.im_dim[ii]);
			vec.rhs_os[timestep][ii] = computeConvolution(vec.rhs_os[timestep][ii], g, inputScalars, w_vec, nRekos, ii);
		}
		vec.rhs_os[timestep][ii].eval();
	}
	outputFP.eval();
	return status;
}

// Required as a separate step in the multi-queue/stream case as the unlocks need to be done after everything is done
inline int finalizeBackwardProjectionAF(AF_im_vectors& vec, const scalarStruct& inputScalars, const Weighting& w_vec, const RecMethods& MethodList, af::array& outputFP, af::array& meanBP,
	const af::array& g, ProjectorClass& proj, const uint32_t timestep, const int ii = 0) {
#ifndef CPU
	const bool fastStep = inputScalars.fastPDHG && proj.fastStep != 0;
#else
	// No CPU support
	const bool fastStep = false;
#endif
#ifndef CPU
	// The finalization is done only once all queues/streams have finished when using multi-resolution
	if (fastStep) {
		if (proj.fastAlg == 0)
			vec.uCP[timestep][ii].unlock();
		vec.im_os[timestep][ii].unlock();
		if (w_vec.precondTypeIm[0] || w_vec.precondTypeIm[1])
			w_vec.D[0][ii].unlock();
#if defined(CUDA) || defined(HIP)
		// The estimate texture is only created for the main volume, so destroy it once
		if (ii == 0 && proj.fastNLMUsed && (inputScalars.FPType == 5 || inputScalars.largeDim)) {
			// The side-stream backprojections that sample this texture may still be in flight; after joinSideQueues the main stream waits for them
			getErrorString(cuStreamSynchronize(proj.CLCommandQueue[0]));
			CUresult statusTex = cuTexObjectDestroy(proj.d_inputI);
			if (statusTex != CUDA_SUCCESS)
				getErrorString(statusTex);
			statusTex = cuArrayDestroy(proj.imArray);
			if (statusTex != CUDA_SUCCESS)
				getErrorString(statusTex);
		}
#endif
	}
	// Gated exactly like the transfer in backwardProjectionAFOpenCL: the patch weights are a shared array,
	// locked once for the main volume and NOT tied to fastStep, so release them under the same condition.
	if (ii == 0 && inputScalars.fastPDHG && proj.fastNLMUsed && MethodList.NLM)
		w_vec.gaussianNLM.unlock();
#endif
	// fastPDHG binds im_os as the backprojection output, so rhs_os was never locked
	if (!fastStep)
		vec.rhs_os[timestep][ii].unlock();
	outputFP.unlock();
	if (inputScalars.meanBP && inputScalars.BPType == 5)
		meanBP.unlock();
	if (!fastStep) {
#if defined(OPENCL) || defined(METAL)
		if (inputScalars.atomic_64bit)
			vec.rhs_os[timestep][ii] = vec.rhs_os[timestep][ii].as(f32) / TH;
		else if (inputScalars.atomic_32bit)
			vec.rhs_os[timestep][ii] = vec.rhs_os[timestep][ii].as(f32) / TH32;
#endif
		if (inputScalars.use_psf) {
			const int nRekos = vec.rhs_os[timestep][ii].elements() / (inputScalars.im_dim[ii]);
			vec.rhs_os[timestep][ii] = computeConvolution(vec.rhs_os[timestep][ii], g, inputScalars, w_vec, nRekos, ii);
		}
		vec.rhs_os[timestep][ii].eval();
	}
	return 0;
}

// Metal does not yet expose the auxiliary projector-side prior kernels
// other than PDHG. Keep their API available so the common reconstruction
// driver reports an unsupported selection cleanly.
#ifdef METAL
inline int MRPAF(af::array&, af::array&, const scalarStruct&, ProjectorClass&, const uint32_t, const uint32_t, const uint32_t) { return -1; }
inline int NLMAF(af::array&, const af::array&, const scalarStruct&, Weighting&, ProjectorClass&, const float, const int = 0) { return -1; }
inline int RDPAF(af::array&, const af::array&, const scalarStruct&, const float, ProjectorClass&, const float, const af::array&, const Weighting&,
	const bool = false, const bool = false, const int = 0) { return -1; }
inline int GGMRFAF(af::array&, const af::array&, const scalarStruct&, const float, const float, const float, const float, ProjectorClass&,
	const float, const Weighting&, const int = 0) { return -1; }
inline int TVAF(af::array&, const af::array&, const scalarStruct&, const float, const TVdata&, ProjectorClass&, const float, const Weighting&, const int = 0) { return -1; }
inline int hyperAF(af::array&, const af::array&, const scalarStruct&, const float, ProjectorClass&, const float, const Weighting&, const int = 0) { return -1; }
inline int proxQAF(std::vector<af::array>&, const float, ProjectorClass&) { return -1; }
inline int proxTVQAF(std::vector<af::array>&, const float, ProjectorClass&) { return -1; }
inline int proxTGVQAF(std::vector<af::array>&, const scalarStruct&, const float, ProjectorClass&) { return -1; }
inline int proxTVDivAF(const std::vector<af::array>&, af::array&, const scalarStruct&, ProjectorClass&, const uint32_t) { return -1; }
inline int proxTVGradAF(const af::array&, std::vector<af::array>&, const scalarStruct&, const float, const std::vector<af::array>&, ProjectorClass&) { return -1; }
inline int proxTGVSymmDerivAF(const std::vector<af::array>&, std::vector<af::array>&, const scalarStruct&, const float, ProjectorClass&) { return -1; }
inline int proxTGVDivAF(const std::vector<af::array>&, std::vector<af::array>&, std::vector<af::array>&, const scalarStruct&,
	const float, const float, ProjectorClass&) { return -1; }
inline int poissonUpdateAF(af::array&, const af::array&, const scalarStruct&, const float, const float, const float, ProjectorClass&, const int = 0) { return -1; }
#else

// Computes custom median root prior, OpenCL, CUDA or CPU
inline int MRPAF(af::array& padd, af::array& grad, const scalarStruct& inputScalars, ProjectorClass& proj, const uint32_t medx, const uint32_t medy, const uint32_t medz) {

	int status = 0;
	if (DEBUG) {
		mexPrintBase("padd = %f\n", af::sum<float>(padd));
	}
	proj.d_W = transferAF(grad);
	proj.d_inputB = transferAF(padd);
#ifndef CPU
	uint64_t global_size[3] = { static_cast<uint64_t>(padd.dims(0)), static_cast<uint64_t>(padd.dims(1)), static_cast<uint64_t>(padd.dims(2)) };
	status = proj.computeMRP(inputScalars, global_size);
	grad.unlock();
	padd.unlock();
	if (status != 0) {
		return -1;
	}
#else
	proj.computeMRP(inputScalars, medx, medy, medz);
	grad.unlock();
	padd.unlock();
#endif
	return status;
}

// NLM
inline int NLMAF(af::array& grad, const af::array& im, const scalarStruct& inputScalars, Weighting& w_vec, ProjectorClass& proj, const float beta, const int kk = 0) {

	int status = 0;
#ifdef CPU
	// Accumulate the NLM gradient to the backprojection
	const af::array dUorig = grad.copy();
	grad = padding(grad, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], w_vec.Ndx + w_vec.Nlx, w_vec.Ndy + w_vec.Nly, w_vec.Ndz + w_vec.Nlz);
	const af::array imApu = padding(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], w_vec.Ndx + w_vec.Nlx, w_vec.Ndy + w_vec.Nly, w_vec.Ndz + w_vec.Nlz);
#endif
	proj.d_W = transferAF(grad);
	proj.d_gaussianNLM = transferAF(w_vec.gaussianNLM);
#ifndef CPU
	if (!inputScalars.useImages) {
		proj.d_inputB = transferAF(im);
	}
	// FPType == 5 caches integral images, not the estimate, so the prior still needs its own copy
	else if (inputScalars.FPType == 5) {
#if defined(CUDA) || defined(HIP)
		uint32_t Nz = inputScalars.Nz[0];
		if (inputScalars.largeDim)
			Nz = inputScalars.lDimStruct.NzPr[kk];
		CUdeviceptr* input = transferAF(im);
		status = proj.transferTex(inputScalars, input, false, Nz);
#else
		if (inputScalars.largeDim) {
			proj.region[2] = inputScalars.lDimStruct.NzPr[kk];
			proj.d_inputI = cl::Image3D(proj.CLContext, CL_MEM_READ_ONLY, proj.format, proj.region[0], proj.region[1], proj.region[2], 0, 0, NULL, &status);
			OCL_CHECK(status, "Failed to create prior image\n", -1);
			if (DEBUG) {
				mexPrintBase("inputScalars.lDimStruct.NzPr[kk] = %u\n", inputScalars.lDimStruct.NzPr[kk]);
				mexPrintBase("inputScalars.lDimStruct.Nz[kk] = %u\n", inputScalars.lDimStruct.Nz[kk]);
				mexPrintBase("proj.region[2] = %u\n", proj.region[2]);
				mexPrintBase("im.elements() = %u\n", im.elements());
				mexPrintBase("grad.elements() = %u\n", grad.elements());
				mexEval();
			}
		}
		status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(cl::Buffer(*im.device<cl_mem>(), true), proj.d_inputI, 0, proj.origin, proj.region);
		if (status != 0) {
			getErrorString(status);
			im.unlock();
			grad.unlock();
			w_vec.gaussianNLM.unlock();
			mexPrint("Failed to copy NLM image\n");
			return -1;
		}
		if (inputScalars.largeDim)
			proj.region[2] = inputScalars.Nz[0];
#endif
	}
	status = proj.computeNLM(inputScalars, w_vec, beta, kk);
	grad.unlock();
	if (!inputScalars.useImages || inputScalars.FPType == 5)
		im.unlock();
	w_vec.gaussianNLM.unlock();
	if (status != 0) {
		return -1;
	}
#else
	proj.d_inputB = transferAF(imApu);
	status = proj.computeNLM(inputScalars, w_vec, beta);
	grad.unlock();
	imApu.unlock();
	w_vec.gaussianNLM.unlock();
	grad = grad(af::seq(w_vec.Ndx + w_vec.Nlx, inputScalars.Nx[0] + w_vec.Ndx + w_vec.Nlx - 1), af::seq(w_vec.Ndy + w_vec.Nly, inputScalars.Ny[0] + w_vec.Ndy + w_vec.Nly - 1), af::seq(w_vec.Ndz + w_vec.Nlz, inputScalars.Nz[0] + w_vec.Ndz + w_vec.Nlz - 1));
	grad = af::flat(grad);
	grad *= beta;
	grad += dUorig;
	grad.eval();
#endif
#ifdef OPENCL
	if (inputScalars.largeDim)
		proj.region[2] = inputScalars.Nz[0];
#endif
	return status;
}

// RDP
inline int RDPAF(af::array& grad, const af::array& im, const scalarStruct& inputScalars, const float gamma, ProjectorClass& proj, const float beta, const af::array& RDPref, const Weighting& w_vec,
	const bool RDPLargeNeighbor = false, const bool useRDPRef = false, const int kk = 0) {
	im.eval();
	int status = 0;
	proj.d_W = transferAF(grad);
#ifndef CPU
	if (inputScalars.useImages) {
#if defined(CUDA) || defined(HIP)
		// FPType == 5 caches integral images, not the estimate, so the prior still needs its own copy
		if (inputScalars.FPType == 5) {
			uint32_t Nz = inputScalars.Nz[0];
			if (inputScalars.largeDim)
				Nz = inputScalars.lDimStruct.NzPr[kk];
			CUdeviceptr* input = transferAF(im);
			status = proj.transferTex(inputScalars, input, false, Nz);
		}
		if (RDPLargeNeighbor && useRDPRef) {
			CUdeviceptr* inputRef = transferAF(RDPref);
			status = proj.transferTex(inputScalars, inputRef, true);
		}
#else
		// The largeDim subvolume must match the predetermined (padded) size
		if (inputScalars.largeDim)
			proj.region[2] = inputScalars.lDimStruct.NzPr[kk];
		// FPType == 5 caches integral images, not the estimate, so the prior still needs its own copy
		if (inputScalars.FPType == 5) {
			if (inputScalars.largeDim) {
				proj.d_inputI = cl::Image3D(proj.CLContext, CL_MEM_READ_ONLY, proj.format, proj.region[0], proj.region[1], proj.region[2], 0, 0, NULL, &status);
				OCL_CHECK(status, "Failed to create prior image\n", -1);
			}
			status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(cl::Buffer(*im.device<cl_mem>(), true), proj.d_inputI, 0, proj.origin, proj.region);
			if (status != 0) {
				getErrorString(status);
				im.unlock();
				grad.unlock();
				mexPrint("Failed to copy RDP image\n");
				return -1;
			}
		}
		if (RDPLargeNeighbor && useRDPRef) {
			status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(cl::Buffer(*RDPref.device<cl_mem>(), true), proj.d_RDPrefI, 0, proj.origin, proj.region);
			if (status != 0) {
				getErrorString(status);
				if (inputScalars.FPType == 5)
					im.unlock();
				grad.unlock();
				mexPrint("Failed to copy RDP image\n");
				return -1;
			}
		}
#endif
	}
	else {
		proj.d_inputB = transferAF(im);
		if (RDPLargeNeighbor && useRDPRef)
			proj.d_RDPref = transferAF(RDPref);
	}
	if (DEBUG) {
		mexPrintBase("im.elements() = %u\n", im.elements());
		mexPrintBase("sum(isnan(im)) = %f\n", af::sum<float>(isNaN(im)));
		mexEval();
	}
	status = proj.computeRDP(inputScalars, gamma, w_vec, beta, kk, RDPLargeNeighbor, useRDPRef);
	grad.unlock();
	if (!inputScalars.useImages || inputScalars.FPType == 5)
		im.unlock();
	if (RDPLargeNeighbor && useRDPRef)
		RDPref.unlock();
	if (status != 0) {
		return -1;
	}
#else
	proj.d_inputB = transferAF(im);
	status = proj.computeRDP(inputScalars, gamma, beta, RDPLargeNeighbor, useRDPRef);
	grad.unlock();
	im.unlock();
#endif
#ifdef OPENCL
	if (inputScalars.largeDim)
		proj.region[2] = inputScalars.Nz[0];
#endif
	return status;
}

// GGMRF
inline int GGMRFAF(af::array& grad, const af::array& im, const scalarStruct& inputScalars, const float p, const float q, const float c, const float pqc, ProjectorClass& proj,
	const float beta, const Weighting& w_vec, const int kk = 0) {
	im.eval();
	int status = 0;
	proj.d_W = transferAF(grad);
#ifndef CPU
	if (!inputScalars.useImages) {
		proj.d_inputB = transferAF(im);
	}
	// FPType == 5 caches integral images, not the estimate, so the prior still needs its own copy
	else if (inputScalars.FPType == 5) {
#if defined(CUDA) || defined(HIP)
		uint32_t Nz = inputScalars.Nz[0];
		if (inputScalars.largeDim)
			Nz = inputScalars.lDimStruct.NzPr[kk];
		CUdeviceptr* input = transferAF(im);
		status = proj.transferTex(inputScalars, input, false, Nz);
#else
		if (inputScalars.largeDim) {
			proj.region[2] = inputScalars.lDimStruct.NzPr[kk];
			proj.d_inputI = cl::Image3D(proj.CLContext, CL_MEM_READ_ONLY, proj.format, proj.region[0], proj.region[1], proj.region[2], 0, 0, NULL, &status);
			OCL_CHECK(status, "Failed to create prior image\n", -1);
		}
		status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(cl::Buffer(*im.device<cl_mem>(), true), proj.d_inputI, 0, proj.origin, proj.region);
		if (status != 0) {
			getErrorString(status);
			im.unlock();
			grad.unlock();
			mexPrint("Failed to copy GGMRF image\n");
			return -1;
		}
#endif
	}
	if (DEBUG) {
		mexPrintBase("im.elements() = %u\n", im.elements());
		mexPrintBase("sum(isnan(im)) = %f\n", af::sum<float>(isNaN(im)));
		mexEval();
	}
	status = proj.computeGGMRF(inputScalars, p, q, c, pqc, w_vec, beta, kk);
	grad.unlock();
	if (!inputScalars.useImages || inputScalars.FPType == 5)
		im.unlock();
	if (status != 0) {
		return -1;
	}
#else
	proj.d_inputB = transferAF(im);
	status = proj.computeGGMRF(inputScalars, p, q, c, pqc, beta);
	grad.unlock();
	im.unlock();
#endif
#ifdef OPENCL
	if (inputScalars.largeDim)
		proj.region[2] = inputScalars.Nz[0];
#endif
	return status;
}

// TV
inline int TVAF(af::array& grad, const af::array& im, const scalarStruct& inputScalars, const float sigma, const TVdata& data, ProjectorClass& proj, const float beta, const Weighting& w_vec, const int kk = 0) {
	im.eval();
	int status = 0;
	int type = 0;
	proj.d_W = transferAF(grad);
	if (data.TV_use_anatomical)
		proj.d_refIm = transferAF(data.refIm);
	float C = 0.f;
	if (data.TVtype == 5) {
		C = data.eta;
		type = 3;
	}
	else if (data.TVtype == 2) {
		C = data.C;
		type = 2;
	}
	else if ((data.TVtype == 1 && data.TV_use_anatomical))
		type = 1;
#ifndef CPU
	if (!inputScalars.useImages) {
		proj.d_inputB = transferAF(im);
	}
	// FPType == 5 caches integral images, not the estimate, so the prior still needs its own copy
	else if (inputScalars.FPType == 5) {
#if defined(CUDA) || defined(HIP)
		uint32_t Nz = inputScalars.Nz[0];
		if (inputScalars.largeDim)
			Nz = inputScalars.lDimStruct.NzPr[kk];
		CUdeviceptr* input = transferAF(im);
		status = proj.transferTex(inputScalars, input, false, Nz);
#else
		if (inputScalars.largeDim) {
			proj.region[2] = inputScalars.lDimStruct.NzPr[kk];
			proj.d_inputI = cl::Image3D(proj.CLContext, CL_MEM_READ_ONLY, proj.format, proj.region[0], proj.region[1], proj.region[2], 0, 0, NULL, &status);
			OCL_CHECK(status, "Failed to create prior image\n", -1);
		}
		status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(cl::Buffer(*im.device<cl_mem>(), true), proj.d_inputI, 0, proj.origin, proj.region);
		if (status != 0) {
			getErrorString(status);
			im.unlock();
			grad.unlock();
			mexPrint("Failed to copy TV image\n");
			return -1;
		}
#endif
	}
	if (DEBUG) {
		mexPrintBase("im.elements() = %u\n", im.elements());
		mexPrintBase("sum(isnan(im)) = %f\n", af::sum<float>(isNaN(im)));
		mexEval();
	}
	const float smooth = data.TVsmoothing;
	status = proj.TVGradient(inputScalars, sigma, smooth, w_vec, beta, kk, C, type);
	grad.unlock();
	if (!inputScalars.useImages || inputScalars.FPType == 5)
		im.unlock();
	if (data.TV_use_anatomical)
		data.refIm.unlock();
	if (status != 0) {
		return -1;
	}
#else
	proj.d_inputB = transferAF(im);
	status = proj.TVGradient(inputScalars, data, sigma, data.TVsmoothing, beta, C, data.TVtype);
	grad.unlock();
	im.unlock();
	if (data.TV_use_anatomical)
		data.refIm.unlock();
#endif
#ifdef OPENCL
	if (inputScalars.largeDim)
		proj.region[2] = inputScalars.Nz[0];
#endif
	return status;
}

// Hyperbolic prior
inline int hyperAF(af::array& grad, const af::array& im, const scalarStruct& inputScalars, const float sigma, ProjectorClass& proj, const float beta, const Weighting& w_vec, const int kk = 0) {
	im.eval();
	int status = 0;
	int type = 0;
#ifndef CPU
	proj.d_W = transferAF(grad);
	if (!inputScalars.useImages) {
		proj.d_inputB = transferAF(im);
	}
	// FPType == 5 caches integral images, not the estimate, so the prior still needs its own copy
	else if (inputScalars.FPType == 5) {
#if defined(CUDA) || defined(HIP)
		uint32_t Nz = inputScalars.Nz[0];
		if (inputScalars.largeDim)
			Nz = inputScalars.lDimStruct.NzPr[kk];
		CUdeviceptr* input = transferAF(im);
		status = proj.transferTex(inputScalars, input, false, Nz);
#else
		if (inputScalars.largeDim) {
			proj.region[2] = inputScalars.lDimStruct.NzPr[kk];
			proj.d_inputI = cl::Image3D(proj.CLContext, CL_MEM_READ_ONLY, proj.format, proj.region[0], proj.region[1], proj.region[2], 0, 0, NULL, &status);
			OCL_CHECK(status, "Failed to create prior image\n", -1);
		}
		status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(cl::Buffer(*im.device<cl_mem>(), true), proj.d_inputI, 0, proj.origin, proj.region);
		if (status != 0) {
			getErrorString(status);
			im.unlock();
			grad.unlock();
			mexPrint("Failed to copy hyperbolic prior image\n");
			return -1;
		}
#endif
	}
	if (DEBUG) {
		mexPrintBase("im.elements() = %u\n", im.elements());
		mexPrintBase("sum(isnan(im)) = %f\n", af::sum<float>(isNaN(im)));
		mexEval();
	}
	status = proj.hyperGradient(inputScalars, sigma, w_vec, beta, kk);
	grad.unlock();
	if (!inputScalars.useImages || inputScalars.FPType == 5)
		im.unlock();
	if (status != 0) {
		return -1;
	}
#else
	proj.d_W = transferAF(grad);
	proj.d_inputB = transferAF(im);
	status = proj.hyperGradient(inputScalars, sigma, beta);
	grad.unlock();
	im.unlock();
#endif
#ifdef OPENCL
	if (inputScalars.largeDim)
		proj.region[2] = inputScalars.Nz[0];
#endif
	return status;
}


#ifndef CPU
// Proximal TV
inline int proxQAF(std::vector<af::array>& q, const float alpha, ProjectorClass& proj) {
	int status = 0;
	uint64_t globalQ = q[0].elements();
	if (DEBUG) {
		mexPrintBase("globalQ = %u\n", globalQ);
		mexPrintBase("q.elements() = %u\n", q[0].elements());
		mexEval();
	}
	proj.d_qX = transferAF(q[0]);
	status = proj.ProxHelperQ(alpha, globalQ);
	q[0].unlock();
	if (status != 0) {
		return -1;
	}
	return status;
}

// Proximal TV
inline int proxTVQAF(std::vector<af::array>& q, const float alpha, ProjectorClass& proj) {
	int status = 0;
	uint64_t globalQ = q[0].elements();
	if (DEBUG) {
		mexPrintBase("globalQ = %u\n", globalQ);
		mexPrintBase("q.elements() = %u\n", q[0].elements());
		mexEval();
	}
	proj.d_qX = transferAF(q[0]);
	proj.d_qY = transferAF(q[1]);
	proj.d_qZ = transferAF(q[2]);
	status = proj.ProxTVHelperQ(alpha, globalQ);
	q[0].unlock();
	q[1].unlock();
	q[2].unlock();
	if (status != 0) {
		return -1;
	}
	return status;
}

// Proximal TGV
inline int proxTGVQAF(std::vector<af::array>& r, const scalarStruct& inputScalars, const float alpha, ProjectorClass& proj) {
	int status = 0;
	uint64_t globalQ = r[0].elements();
	if (DEBUG) {
		mexPrintBase("globalQ = %u\n", globalQ);
		mexPrintBase("q.elements() = %u\n", r[0].elements());
		mexEval();
	}
	proj.d_rX = transferAF(r[0]);
	proj.d_rY = transferAF(r[1]);
	if (!inputScalars.TGV2D) {
		proj.d_rZ = transferAF(r[2]);
		proj.d_rXY = transferAF(r[3]);
		proj.d_rXZ = transferAF(r[4]);
		proj.d_rYZ = transferAF(r[5]);
	}
	else {
		proj.d_rXY = transferAF(r[2]);
	}
	status = proj.ProxTGVHelperQ(inputScalars, alpha, globalQ);
	r[0].unlock();
	r[1].unlock();
	r[2].unlock();
	if (!inputScalars.TGV2D) {
		r[3].unlock();
		r[4].unlock();
		r[5].unlock();
	}
	if (status != 0) {
		return -1;
	}
	return status;
}

// Proximal TV divergence
inline int proxTVDivAF(const std::vector<af::array>& grad, af::array& input, const scalarStruct& inputScalars, ProjectorClass& proj, const uint32_t timestep) {
	int status = 0;
	if (DEBUG) {
		mexPrintBase("input.dims(0) = %u\n", input.dims(0));
		mexPrintBase("grad[0].dims(0) = %u\n", grad[0].dims(0));
		mexPrintBase("grad[0].dims(1) = %u\n", grad[0].dims(1));
		mexPrintBase("grad[0].dims(2) = %u\n", grad[0].dims(2));
		mexEval();
	}
	proj.d_qX = transferAF(grad[0]);
	proj.d_qY = transferAF(grad[1]);
	proj.d_qZ = transferAF(grad[2]);
	proj.vec_opencl.d_rhs_os[0] = transferAF(input);
	status = proj.ProxTVDiv(inputScalars);
	grad[0].unlock();
	grad[1].unlock();
	grad[2].unlock();
	input.unlock();
	if (status != 0) {
		return -1;
	}
	return status;
}

// Proximal TV gradient
inline int proxTVGradAF(const af::array& im, std::vector<af::array>& grad, const scalarStruct& inputScalars, const float sigma2, const std::vector<af::array>& v, ProjectorClass& proj) {
	int status = 0;
	if (DEBUG) {
		mexPrintBase("output.dims(0) = %u\n", grad[0].dims(0));
		mexPrintBase("im.dims(0) = %u\n", im.dims(0));
		mexPrintBase("im.dims(1) = %u\n", im.dims(1));
		mexPrintBase("im.dims(2) = %u\n", im.dims(2));
		mexEval();
	}
	const size_t vSize = v.size();
	proj.d_qX = transferAF(grad[0]);
	proj.d_qY = transferAF(grad[1]);
	proj.d_qZ = transferAF(grad[2]);
	proj.d_inputB = transferAF(im);
	if (v.size() > 0) {
		if (DEBUG) {
			mexPrintBase("v0.dims(0) = %u\n", v[0].dims(0));
			mexEval();
		}
		if (DEBUG) {
			mexPrintBase("v1.dims(0) = %u\n", v[1].dims(0));
			mexEval();
		}
		proj.d_vX = transferAF(v[0]);
		proj.d_vY = transferAF(v[1]);
		if (!inputScalars.TGV2D)
			proj.d_vZ = transferAF(v[2]);
	}
	status = proj.ProxTVGrad(inputScalars, sigma2, vSize);
	grad[0].unlock();
	grad[1].unlock();
	grad[2].unlock();
	im.unlock();
	if (v.size() > 0) {
		v[0].unlock();
		v[1].unlock();
		if (!inputScalars.TGV2D) {
			v[2].unlock();
		}
	}
	if (status != 0) {
		return -1;
	}
	return status;
}

// Proximal TGV symmetric derivative
inline int proxTGVSymmDerivAF(const std::vector<af::array>& v, std::vector<af::array>& q, const scalarStruct& inputScalars, const float sigma2, ProjectorClass& proj) {
	int status = 0;
	if (DEBUG) {
		mexPrintBase("input.dims(0) = %u\n", v[0].dims(0));
		if (!inputScalars.TGV2D)
			mexPrintBase("input2.dims(0) = %u\n", v[2].dims(0));
		mexPrintBase("im.dims(0) = %u\n", q[0].dims(0));
		mexPrintBase("im.dims(1) = %u\n", q[0].dims(1));
		mexPrintBase("im.dims(2) = %u\n", q[0].dims(2));
		mexPrintBase("q1.dims(0) = %u\n", q[1].dims(0));
		mexPrintBase("q2.dims(0) = %u\n", q[2].dims(0));
		if (!inputScalars.TGV2D) {
			mexPrintBase("q3.dims(0) = %u\n", q[3].dims(0));
			mexPrintBase("q5.dims(0) = %u\n", q[5].dims(0));
		}
		mexPrintBase("v.size() = %u\n", v.size());
		mexPrintBase("q.size() = %u\n", q.size());
		mexEval();
	}
	proj.d_rX = transferAF(q[0]);
	proj.d_rY = transferAF(q[1]);
	if (!inputScalars.TGV2D) {
		proj.d_rZ = transferAF(q[2]);
		proj.d_rXY = transferAF(q[3]);
		proj.d_rXZ = transferAF(q[4]);
		proj.d_rYZ = transferAF(q[5]);
	}
	else
		proj.d_rXY = transferAF(q[2]);
	proj.d_vX = transferAF(v[0]);
	proj.d_vY = transferAF(v[1]);
	if (!inputScalars.TGV2D)
		proj.d_vZ = transferAF(v[2]);
	status = proj.ProxTGVSymmDeriv(inputScalars, sigma2);
	v[0].unlock();
	v[1].unlock();
	if (!inputScalars.TGV2D)
		v[2].unlock();
	q[0].unlock();
	q[1].unlock();
	q[2].unlock();
	if (!inputScalars.TGV2D) {
		q[3].unlock();
		q[4].unlock();
		q[5].unlock();
	}
	if (status != 0) {
		return -1;
	}
	return status;
}

// Proximal TGV divergence
inline int proxTGVDivAF(const std::vector<af::array>& q, std::vector<af::array>& v, std::vector<af::array>& p, const scalarStruct& inputScalars,
	const float theta, const float tau, ProjectorClass& proj) {
	int status = 0;
	if (DEBUG) {
		mexPrintBase("v.dims(0) = %u\n", v[0].dims(0));
		mexPrintBase("q2.dims(0) = %u\n", q[0].dims(0));
		mexPrintBase("q2.dims(1) = %u\n", q[0].dims(1));
		mexPrintBase("q2.dims(2) = %u\n", q[0].dims(2));
		mexEval();
	}
	proj.d_rX = transferAF(q[0]);
	proj.d_rY = transferAF(q[1]);
	if (!inputScalars.TGV2D) {
		proj.d_rZ = transferAF(q[2]);
		proj.d_rXY = transferAF(q[3]);
		proj.d_rXZ = transferAF(q[4]);
		proj.d_rYZ = transferAF(q[5]);
	}
	else
		proj.d_rXY = transferAF(q[2]);
	proj.d_vX = transferAF(v[0]);
	proj.d_vY = transferAF(v[1]);
	if (!inputScalars.TGV2D)
		proj.d_vZ = transferAF(v[2]);
	proj.d_qX = transferAF(p[0]);
	proj.d_qY = transferAF(p[1]);
	proj.d_qZ = transferAF(p[2]);
	status = proj.ProxTGVDiv(inputScalars, theta, tau);
	v[0].unlock();
	v[1].unlock();
	if (!inputScalars.TGV2D)
		v[2].unlock();
	q[0].unlock();
	q[1].unlock();
	q[2].unlock();
	if (!inputScalars.TGV2D) {
		q[3].unlock();
		q[4].unlock();
		q[5].unlock();
	}
	p[0].unlock();
	p[1].unlock();
	p[2].unlock();
	if (status != 0) {
		return -1;
	}
	return status;
}

// Poisson-estimate update kernel call for PKMA, MBSREM and BSREM
inline int poissonUpdateAF(af::array& im, const af::array& rhs, const scalarStruct& inputScalars, const float lambda, const float epps, const float alpha, ProjectorClass& proj, const int ii = 0) {
	int status = 0;
	proj.d_im = transferAF(im);
	proj.d_rhs = transferAF(rhs);
	status = proj.PoissonUpdate(inputScalars, lambda, epps, alpha, ii);
	rhs.unlock();
	im.unlock();
	if (status != 0) {
		return -1;
	}
	return status;
}

#endif
#endif

#ifndef CPU
// Computes elementwise multiplication or division, depending on mult-variable
inline int elementWiseAF(const af::array& vector, af::array& input, const bool mult, ProjectorClass& proj, const bool D2 = false) {
	int status = 0;
	uint64_t gSize[3] = { static_cast<uint64_t>(input.dims(0)), static_cast<uint64_t>(input.dims(1)), static_cast<uint64_t>(input.dims(2)) };
	if (DEBUG) {
		mexPrintBase("vector.dims[0] = %u\n", vector.dims(0));
		mexPrintBase("vector.dims[1] = %u\n", vector.dims(1));
		mexPrintBase("vector.dims[2] = %u\n", vector.dims(2));
		mexPrintBase("input.dims[0] = %u\n", input.dims(0));
		mexPrintBase("input.dims[1] = %u\n", input.dims(1));
		mexPrintBase("input.dims[2] = %u\n", input.dims(2));
		mexEval();
	}
	proj.d_vector = transferAF(vector);
	proj.d_input = transferAF(input);
	status = proj.elementWiseComp(mult, gSize, D2);
	vector.unlock();
	input.unlock();
	if (status != 0) {
		return -1;
	}
	return status;
}

// Same as above, but for PDHG
inline int PDHGUpdateAF(af::array& im, const af::array& rhs, const scalarStruct& inputScalars, AF_im_vectors& vec, const float epps, const float theta, const float tau, ProjectorClass& proj, const uint32_t timestep, const int ii = 0) {
	int status = 0;
	proj.d_im = transferAF(im);
	proj.d_rhs = transferAF(rhs);
	proj.d_U = transferAF(vec.uCP[timestep][ii]);
	status = proj.PDHGUpdate(inputScalars, epps, theta, tau, ii);
	rhs.unlock();
	im.unlock();
	vec.uCP[timestep][ii].unlock();
	if (status != 0) {
		return -1;
	}
	return status;
}
#endif

#if !defined(CPU)
#if defined(METAL)
inline int rotateCustomAF(af::array& imrot, const af::array& im, const scalarStruct& inputScalars, ProjectorClass& proj, const float angle, const int ii = 0) {
	int status = 0;
	if (!inputScalars.useBuffers) {
#if defined(CUDA) || defined(HIP) || defined(METAL)
		AFDEVBUFF_t input = transferAF(im);
		status = proj.transferTex(inputScalars, input, false, inputScalars.Nz[0]);
#elif defined(OPENCL)
        status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(cl::Buffer(*im.device<cl_mem>(), true), proj.d_inputI, 0, proj.origin, proj.region);
        if (status != 0) {
			getErrorString(status);
			im.unlock();
			mexPrint("Failed to copy rotation image\n");
			return -1;
		}
#endif
    } else {
        proj.d_im = transferAF(im);
    }
	proj.d_rhs = transferAF(imrot);
	const float cosa = std::cos(-angle);
	const float sina = std::sin(-angle);
	status = proj.rotateCustom(inputScalars, cosa, sina, ii);
	imrot.unlock();
	im.unlock();
	if (status != 0) {
		return -1;
	}
	return status;
}
#else
// Rotates the input image and, at the same time, shifts it by (offX, offY) voxels
// The output is obtained as imrot(q) = im(R * q + off), where q is the (centered) output coordinate
// and R the rotation matrix. Combining the rotation and the translation means that the image is
// interpolated only once, which matters since the transaxial CDRF is compensated for the blurring
// caused by a single interpolation
// texSlot selects which of the projector's rotation images is used (0 = emission volume, 1 = attenuation
// volume) and upload can be set to false when the image is already there, which is the case for the
// forward projection where the same image is rotated once for every projection, and for the
// attenuation map which never changes
// With accumulate the rotated image is added to imrot, which then has to hold the output already (only its
// device pointer is used), and the elements that the rotation does not reach are left untouched
inline int rotateCustomAF(af::array& imrot, const af::array& im, const scalarStruct& inputScalars, ProjectorClass& proj, const float angle, const float offX, const float offY,
	RotWindow win, const int ii = 0, const int texSlot = 0, const bool upload = true, const bool accumulate = false) {
	int status = 0;
	// Without a window in x (the depth) the source and the output are whole in that dimension
	if (win.outNxW <= 0)
		win.outNxW = win.outNx;
	if (win.srcNxB <= 0)
		win.srcNxB = win.srcNx;
	if (!inputScalars.useBuffers) {
		if (upload) {
			// Only the part of the image that holds data is uploaded, the rest of the image is not read
			// HIP does not define CUDA, so it has to be listed separately or the image is never uploaded
#if defined(CUDA) || defined(HIP)
			AFDEVBUFF_t input = transferAF(im);
			status = proj.transferRotTex(inputScalars, input, win, texSlot);
			im.unlock();
			if (status != 0)
				return -1;
#elif defined(OPENCL)
			if (proj.ensureRotImage(win, texSlot) != 0)
				return -1;
			const cl::detail::size_t_array rotOrigin = { { static_cast<::size_t>(win.srcX), static_cast<::size_t>(win.srcY), static_cast<::size_t>(win.srcZ) } };
			const cl::detail::size_t_array rotRegion = { { static_cast<::size_t>(win.srcNxB), static_cast<::size_t>(win.srcNy), static_cast<::size_t>(win.srcNz) } };
			status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(cl::Buffer(*im.device<cl_mem>(), true), proj.d_rotI[texSlot], 0, rotOrigin, rotRegion);
			im.unlock();
			if (status != 0) {
				getErrorString(status);
				mexPrint("Failed to copy rotation image\n");
				return -1;
			}
#endif
		}
    } else {
        proj.d_im = transferAF(im);
    }
	proj.d_rhs = transferAF(imrot);
	const float cosa = std::cos(-angle);
	const float sina = std::sin(-angle);
	status = proj.rotateCustom(inputScalars, cosa, sina, offX, offY, win, texSlot, accumulate ? 1 : 0);
	imrot.unlock();
	if (inputScalars.useBuffers)
		im.unlock();
	if (status != 0) {
		return -1;
	}
	return status;
}

// One pass of the separable collimator-detector response convolution. The filter bank holds one 1D
// filter for each slice of the input along its third dimension
inline int convolveCDRFAF(af::array& output, const af::array& input, const af::array& filter, const scalarStruct& inputScalars, ProjectorClass& proj, const int axis) {
	proj.d_im = transferAF(input);
	proj.d_rhs = transferAF(output);
	proj.d_convFilter = transferAF(filter);
	const int fRadius = static_cast<int>((filter.dims(0) * filter.dims(1) - 1) / 2);
	const int status = proj.convolveCDRF(inputScalars, static_cast<int>(input.dims(0)), static_cast<int>(input.dims(1)), static_cast<int>(input.dims(2)), fRadius, axis);
	input.unlock();
	output.unlock();
	filter.unlock();
	return status;
}
#endif // END METAL
#endif

// Various batch functions
inline af::array batchMinus(const af::array& lhs, const af::array& rhs) {
	return lhs - rhs;
}

inline af::array batchPlus(const af::array& lhs, const af::array& rhs) {
	return lhs + rhs;
}

inline af::array batchMul(const af::array& lhs, const af::array& rhs) {
	return lhs * rhs;
}

inline af::array batchDiv(const af::array& lhs, const af::array& rhs) {
	return lhs / rhs;
}

inline af::array batchNotEqual(const af::array& lhs, const af::array& rhs) {
	return lhs != rhs;
}

// Compute the epsilon value for the MBSREM/MRAMLA
inline float MBSREM_epsilon(const af::array& Sino, const af::array& D, const float epps = 1e-5f, const uint32_t randoms_correction = false, const af::array& rand = af::constant(0.1, 1, 1),
	const bool TOF = false, const int64_t nBins = 1, const bool CT = false) {
	float eps;
	if (CT) {
		af::array hk_summa = -af::exp(-Sino) / Sino - Sino;
		hk_summa(af::isNaN(hk_summa)) = 0.f;
		af::array P_Sino, apu, Iind;
		if (randoms_correction == 1u) {
			Iind = (Sino > 0.f & rand == 0.f);
			if (af::sum<float>(Iind) == 0.f)
				return 1e8f;
			P_Sino = Sino(Iind);
			apu = D + rand;
			apu = af::sum(-af::exp(-apu) / Sino - apu);
		}
		else {
			Iind = (Sino > 0.f);
			P_Sino = Sino(Iind);
			apu = af::sum(-af::exp(-D) / Sino - D);
		}
		if (randoms_correction == 1u)
			hk_summa = af::batchFunc(af::sum(hk_summa), hk_summa(Iind), batchMinus);
		else
			hk_summa = af::batchFunc(af::sum(hk_summa), hk_summa(Iind), batchMinus);
		af::array epsilon = (af::min)(P_Sino, af::log(af::batchFunc(apu, hk_summa, batchMinus) / P_Sino));
		eps = af::min<float>(epsilon);
	}
	else {
		af::array hk_summa = Sino * af::log(Sino) - Sino;
		hk_summa(af::isNaN(hk_summa)) = 0.f;
		af::array P_Sino, apu, Iind;
		if (TOF && randoms_correction) {
			af::array rInd = rand == 0.f;
			P_Sino = Sino(Sino > 0.f & af::tile(rInd, nBins));
			apu = D + af::tile(rand, nBins);
			apu = af::sum(Sino * af::log(apu) - apu);
			hk_summa = af::batchFunc(af::sum(hk_summa), hk_summa(Sino > 0.f & af::tile(rInd, nBins)), batchMinus);
		}
		else {
			if (randoms_correction == 1u) {
				Iind = (Sino > 0.f & rand == 0.f);
				if (af::sum<float>(Iind) == 0.f)
					return 1e8f;
				P_Sino = Sino(Iind);
				apu = D + rand;
				apu = af::sum(Sino * af::log(apu) - apu);
			}
			else {
				Iind = (Sino > 0.f);
				P_Sino = Sino(Iind);
				apu = af::sum(Sino * af::log(D) - D);
			}
			if (randoms_correction == 1u)
				hk_summa = af::batchFunc(af::sum(hk_summa), hk_summa(Iind), batchMinus);
			else
				hk_summa = af::batchFunc(af::sum(hk_summa), hk_summa(Iind), batchMinus);
		}
		af::array epsilon = (af::min)(P_Sino, af::exp(af::batchFunc(apu, hk_summa, batchMinus) / P_Sino));
		eps = af::min<float>(epsilon);
	}
	eps = eps <= 0.f ? epps : eps;
	return eps;
}

// An array that owns its (contiguous) buffer, for the sources that every view of a projector call hands to a kernel as a device pointer
// (the emission volumes, the measurements). ArrayFire's device() makes a copy of the whole array every time it is called on a handle
// that shares its buffer with another one (af::moddims of a stored image, an array that is also kept in a vector or by the caller), is
// a sub-array or does not own its buffer, and that copy is lost again at the unlock. Done once here, the pointer of the array is
// the same for every view. Nothing is done without a device (CPU), where there is no device pointer
inline af::array ownedType6(const af::array& a) {
#ifdef CPU
	return a;
#else
	af::array c = a.copy();
	c.eval();
	return c;
#endif
}

// Rounds n up to a multiple of q (n if q <= 1). The arrays that change size from one view to the next are allocated with their
// dimensions rounded up, so that the memory pool of ArrayFire finds the buffer of the previous view again instead of allocating a
// new one and freeing the old one, which is what the sizes that follow the window of the view would do (see RotWindow::allocQ). The
// kernels get the logical sizes and the allocated dimensions as strides, so the elements that the rounding adds are never read
inline int roundUpType6(const int n, const int q) {
	return q > 1 ? ((n + q - 1) / q) * q : n;
}
// The multiple that the arrays of the projector that works on every volume of a view together are rounded up to
static const int type6AllocQuantum = 16;

#if defined(METAL)
inline af::array rotateHelperType6(const af::array& im, const scalarStruct& inputScalars, ProjectorClass& proj, const float angleDeg, const int ii = 0) {
    af::array imrot;
    const float angleRad = angleDeg * static_cast<float>(M_PI) / 180.f;
#ifndef CPU
    imrot = af::constant(0.f, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii]);
    rotateCustomAF(imrot, im, inputScalars, proj, angleRad, ii);
#else
    imrot = af::rotate(im, angleRad, true, AF_INTERP_BILINEAR);
#endif
    return imrot;
}
#else
// Rotates the image by angleDeg degrees and simultaneously shifts it such that
// imrot(q) = im(R * q + (offX, offY)), where R is the rotation matrix corresponding to angleDeg
// and q the (centered) output coordinate. The offsets are given in voxels and may be non-integer,
// which is required with swivelling detector heads where the detector panel is, in general, not
// centered with respect to the rotation axis
// win selects which part of the volume is written and which part of the input holds data. Both default
// to the whole volume, which gives a plain rotation
inline af::array rotateHelperType6(const af::array& im, const scalarStruct& inputScalars, ProjectorClass& proj, const float angleDeg, const float offX = 0.f, const float offY = 0.f, const int ii = 0,
    const int texSlot = 0, const bool upload = true, RotWindow win = RotWindow()) {
    af::array imrot;
    const float angleRad = angleDeg * static_cast<float>(M_PI) / 180.f;
    // Without a window the input volume is also the output volume, which is a plain rotation
    if (win.srcNx <= 0) {
        win.srcNx = static_cast<int>(inputScalars.Nx[ii]);
        win.srcNyFull = static_cast<int>(inputScalars.Ny[ii]);
        win.srcNzFull = static_cast<int>(inputScalars.Nz[ii]);
    }
    if (win.outNx <= 0) {
        win.outNx = win.srcNx;
        win.outNyFull = win.srcNyFull;
    }
    // The window in x (the depth), whole unless given
    if (win.outNxW <= 0)
        win.outNxW = win.outNx;
    if (win.srcNxB <= 0)
        win.srcNxB = win.srcNx;
    if (win.outNy <= 0) {
        win.outNy = win.outNyFull;
        win.outNz = win.srcNzFull;
    }
    if (win.srcNy <= 0) {
        win.srcNy = win.srcNyFull;
        win.srcNz = win.srcNzFull;
    }
#ifndef CPU
    // The kernel writes every element of the window of the output, so it does not have to be zeroed here. With allocQ the array is
    // allocated larger than that (the elements outside of the window are not written, nor read by whoever uses the array)
    if (win.allocQ > 1) {
        win.outSx = roundUpType6(win.outNxW, win.allocQ);
        win.outSy = roundUpType6(win.outNy, win.allocQ);
        imrot = af::array(af::dim4(win.outSx, win.outSy, roundUpType6(win.outNz, win.allocQ)), f32);
    }
    else
        imrot = af::array(af::dim4(win.outNxW, win.outNy, win.outNz), f32);
    rotateCustomAF(imrot, im, inputScalars, proj, angleRad, offX, offY, win, ii, texSlot, upload);
#else
    // Translating first and rotating afterwards is equivalent to im(R * q + off). Note that
    // af::translate shifts the sampling point, not the content, so the offsets are used as they are.
    // Unlike the kernel above this interpolates twice, which blurs the result slightly more than the
    // transaxial CDRF assumes, and the two differ in the outermost ring where one of them samples
    // outside the image
    // ArrayFire always works on the whole image, so a banded source is placed back into a full volume
    // and a banded output is cut out afterwards. ArrayFire also rotates the volume into itself, so the
    // multi-resolution volumes, which are rotated into a grid of their own, are not supported here
    af::array input = im;
    if (win.srcNy != win.srcNyFull || win.srcNz != win.srcNzFull) {
        input = af::constant(0.f, win.srcNx, win.srcNyFull, win.srcNzFull);
        input(af::span, af::seq(win.srcY, win.srcY + win.srcNy - 1), af::seq(win.srcZ, win.srcZ + win.srcNz - 1)) = im;
    }
    if (offX != 0.f || offY != 0.f)
        imrot = af::rotate(af::translate(input, offX, offY, 0, 0, AF_INTERP_BILINEAR), angleRad, true, AF_INTERP_BILINEAR);
    else
        imrot = af::rotate(input, angleRad, true, AF_INTERP_BILINEAR);
    if (win.outNy != win.outNyFull || win.outNz != win.srcNzFull)
        imrot = imrot(af::span, af::seq(win.outY, win.outY + win.outNy - 1), af::seq(win.outZ, win.outZ + win.outNz - 1));
#endif
    return imrot;
}

#ifndef CPU
// Rotates and shifts the image like rotateHelperType6, but adds the result to output, which holds the
// whole output window (it is the volume that the backprojection accumulates into), instead of returning
// a new array. What the rotation does not reach is not touched, and nothing else may refer to the
// buffer of output while this runs. The source is uploaded as it is by rotateHelperType6 for texSlot 0
inline void rotateAccumulateType6(af::array& output, const af::array& im, const scalarStruct& inputScalars, ProjectorClass& proj, const float angleDeg, const float offX,
    const float offY, const int ii, const RotWindow& win) {
    const float angleRad = angleDeg * static_cast<float>(M_PI) / 180.f;
    if (static_cast<int64_t>(output.elements()) != static_cast<int64_t>(win.outNxW > 0 ? win.outNxW : win.outNx) * win.outNy * win.outNz) {
        mexPrint("The rotation accumulates into an array of the wrong size");
        return;
    }
    rotateCustomAF(output, im, inputScalars, proj, angleRad, offX, offY, win, ii, 0, true, true);
}
#endif
#endif

inline af::array translateHelperType6(const af::array im, const int shift0, const int shift1, const int shift2) { // TODO: custom translation kernel
    af::array imtrans;
    imtrans = af::shift(im, shift0, shift1, shift2); // ArrayFire performs circular shift.
    if (shift0 > 0) {
        imtrans(af::seq(0, shift0 - 1), af::span, af::span) = 0.f;
    } else if (shift0 < 0) {
        imtrans(af::seq(af::end + shift0 + 1, af::end), af::span, af::span) = 0.f;
    }
    if (shift1 > 0) {
        imtrans(af::span, af::seq(0, shift1 - 1), af::span) = 0.f;
    } else if (shift1 < 0) {
        imtrans(af::span, af::seq(af::end + shift1 + 1, af::end), af::span) = 0.f;
    }
    if (shift2 > 0) {
        imtrans(af::span, af::span, af::seq(0, shift2 - 1)) = 0.f;
    } else if (shift2 < 0) {
        imtrans(af::span, af::span, af::seq(af::end + shift2 + 1, af::end)) = 0.f;
    }
    return imtrans;
}

// Resamples a dense attenuation map (covering the union bounding box of every multi-resolution volume, at volume 0's own, finest, voxel
// size) onto an arbitrary target lattice, using the PHYSICAL position (bx/by/bz), voxel size (dx/dy/dz) and dimensions (nx/ny/nz) of both
//
// Implemented as a mean-preserving box average via a 3D summed-area table (SAT): for each target voxel, the exact physical box it occupies
// is mapped into fractional source-map index space, the (zero-padded) SAT is trilinearly sampled at the box's 8 corners, and the standard
// inclusion-exclusion box-sum is divided by the box's volume (in source-voxel units) to get the mean mu. This is preferred over plain
// trilinear point-sampling of the map because multi-resolution side volumes are normally coarser than the map (multiResolutionScale < 1),
// so a single point sample would alias/under-sample the many source voxels each target voxel actually covers; the box average reduces to
// (numerically) the same result as point-sampling when a target voxel is no coarser than the source (its box collapses towards a single
// corner), so one formula covers both directions.
//
// The table, its trilinear sampling at the box corners and the inclusion-exclusion are all done on the host in double and only the result
// is rounded to float: the prefix sums are ~1e4-1e5 against a voxel of ~0.015, so in float32 the box averages would lose ~1e-3 relative
// accuracy (catastrophic cancellation) and, as the scan order of ArrayFire differs between the backends, the result would differ between
// them. This one is the same on every backend (it needs no f64 support of the device, and it is computed once, as the attenuation source
// is cached for the whole reconstruction). The map atten is NxMap x NyMap x NzMap, the x index the fastest; the result is nx x ny x nz,
// the same order. The box corners are clamped to the map, so the part of a target voxel that is outside the map contributes nothing, i.e.
// no attenuating material; with clampedVolu the box sum is divided by the clamped box (at least 1e-6, the Metal/NumPy convention),
// otherwise by the whole, unclamped, nominal box (which matters for the axis-centered lattice of the rotation-based projector, which
// reaches past the map in the corners)
inline std::vector<float> resampleAttenuationType6Host(const float* atten, const float bxMap, const float byMap, const float bzMap,
	const float dxMap, const float dyMap, const float dzMap, const int NxMap, const int NyMap, const int NzMap,
	const float bxV, const float byV, const float bzV, const float dx, const float dy, const float dz, const int nx, const int ny, const int nz,
	const bool clampedVolu) {
	const size_t sx = static_cast<size_t>(NxMap) + 1U, sy = static_cast<size_t>(NyMap) + 1U, sz = static_cast<size_t>(NzMap) + 1U;
	// Zero-padded 3D summed-area table: sat(i,j,k) = sum of every source voxel with index < (i,j,k)
	std::vector<double> sat(sx * sy * sz, 0.0);
	for (int k = 0; k < NzMap; k++)
		for (int j = 0; j < NyMap; j++) {
			const float* src = atten + (static_cast<size_t>(k) * static_cast<size_t>(NyMap) + static_cast<size_t>(j)) * static_cast<size_t>(NxMap);
			double* dst = sat.data() + ((static_cast<size_t>(k) + 1U) * sy + static_cast<size_t>(j) + 1U) * sx + 1U;
			double run = 0.0;
			for (int i = 0; i < NxMap; i++) {
				run += static_cast<double>(src[i]);
				dst[i] = run;
			}
		}
	for (size_t k = 1; k < sz; k++)
		for (size_t j = 2; j < sy; j++) {
			double* cur = sat.data() + (k * sy + j) * sx;
			const double* prv = cur - sx;
			for (size_t i = 1; i < sx; i++)
				cur[i] += prv[i];
		}
	for (size_t k = 2; k < sz; k++)
		for (size_t j = 0; j < sy; j++) {
			double* cur = sat.data() + (k * sy + j) * sx;
			const double* prv = cur - sx * sy;
			for (size_t i = 0; i < sx; i++)
				cur[i] += prv[i];
		}

	// Physical voxel-boundary edges of the target along each axis, converted into fractional map-index space (0..NMap), clamped to the
	// map, and the lower and the upper source index and the weight of the upper one
	struct Axis { std::vector<double> e, f; std::vector<size_t> i0, i1; };
	auto makeAxis = [](const double b, const double d, const double bMap, const double dMap, const int n, const int nMap) {
		Axis a;
		a.e.resize(static_cast<size_t>(n) + 1U);
		a.f.resize(static_cast<size_t>(n) + 1U);
		a.i0.resize(static_cast<size_t>(n) + 1U);
		a.i1.resize(static_cast<size_t>(n) + 1U);
		for (int i = 0; i <= n; i++) {
			const double e = (std::min)((std::max)((b + static_cast<double>(i) * d - bMap) / dMap, 0.0), static_cast<double>(nMap));
			const double fl = (std::min)((std::max)(std::floor(e), 0.0), static_cast<double>(nMap));
			a.e[i] = e;
			a.i0[i] = static_cast<size_t>(fl);
			a.i1[i] = (std::min)(a.i0[i] + 1U, static_cast<size_t>(nMap));
			a.f[i] = e - fl;
		}
		return a;
	};
	const Axis ax = makeAxis(bxV, dx, bxMap, dxMap, nx, NxMap);
	const Axis ay = makeAxis(byV, dy, byMap, dyMap, ny, NyMap);
	const Axis az = makeAxis(bzV, dz, bzMap, dzMap, nz, NzMap);

	// The summed-area table, trilinearly sampled at the corners (i,j) of the target voxels in the z plane k
	const size_t px = static_cast<size_t>(nx) + 1U, py = static_cast<size_t>(ny) + 1U;
	auto cornerPlane = [&](const size_t k, std::vector<double>& plane) {
		const size_t z0 = az.i0[k], z1 = az.i1[k];
		const double fz = az.f[k];
		for (size_t j = 0; j < py; j++) {
			const size_t y0 = ay.i0[j], y1 = ay.i1[j];
			const double fy = ay.f[j];
			const double* r000 = sat.data() + (z0 * sy + y0) * sx;
			const double* r010 = sat.data() + (z0 * sy + y1) * sx;
			const double* r001 = sat.data() + (z1 * sy + y0) * sx;
			const double* r011 = sat.data() + (z1 * sy + y1) * sx;
			double* out = plane.data() + j * px;
			for (size_t i = 0; i < px; i++) {
				const size_t x0 = ax.i0[i], x1 = ax.i1[i];
				const double fx = ax.f[i];
				const double c00 = r000[x0] * (1.0 - fx) + r000[x1] * fx;
				const double c10 = r010[x0] * (1.0 - fx) + r010[x1] * fx;
				const double c01 = r001[x0] * (1.0 - fx) + r001[x1] * fx;
				const double c11 = r011[x0] * (1.0 - fx) + r011[x1] * fx;
				const double c0 = c00 * (1.0 - fy) + c10 * fy;
				const double c1 = c01 * (1.0 - fy) + c11 * fy;
				out[i] = c0 * (1.0 - fz) + c1 * fz;
			}
		}
	};

	// Nominal size of a target voxel in map voxels. It is the same for every voxel and equal to the clamped size wherever the voxel is
	// inside the map
	const double voluNom = (static_cast<double>(dx) / static_cast<double>(dxMap)) * (static_cast<double>(dy) / static_cast<double>(dyMap))
		* (static_cast<double>(dz) / static_cast<double>(dzMap));
	std::vector<float> res(static_cast<size_t>(nx) * static_cast<size_t>(ny) * static_cast<size_t>(nz));
	std::vector<double> lo(px * py), hi(px * py);
	cornerPlane(0, lo);
	for (int k = 0; k < nz; k++) {
		cornerPlane(static_cast<size_t>(k) + 1U, hi);
		const double dZe = az.e[static_cast<size_t>(k) + 1U] - az.e[static_cast<size_t>(k)];
		for (int j = 0; j < ny; j++) {
			const double dYe = ay.e[static_cast<size_t>(j) + 1U] - ay.e[static_cast<size_t>(j)];
			const size_t r0 = static_cast<size_t>(j) * px, r1 = r0 + px;
			float* out = res.data() + (static_cast<size_t>(k) * static_cast<size_t>(ny) + static_cast<size_t>(j)) * static_cast<size_t>(nx);
			for (int i = 0; i < nx; i++) {
				const size_t a = static_cast<size_t>(i), b = a + 1U;
				const double boxSum =
					  hi[r1 + b] - hi[r1 + a] - hi[r0 + b] - lo[r1 + b]
					+ hi[r0 + a] + lo[r1 + a] + lo[r0 + b] - lo[r0 + a];
				double volu = voluNom;
				if (clampedVolu)
					volu = (std::max)((ax.e[b] - ax.e[a]) * dYe * dZe, 1e-6);
				out[i] = static_cast<float>(boxSum / volu);
			}
		}
		lo.swap(hi);
	}
	return res;
}

// The views of the type-6 projector. The projector is called for the views uu ... uu + length - 1 of a subset, where uu is the running view counter of the
// subsets of ONE frame (it restarts for every frame). The per-view arrays that hold the views of all the frames one frame after the other (the geometry swivelAngles
// and blurPlanes*, the DetectorVector, the per-projection normalization of dynamic data) are read at the offset of the frame in addition, the ones that are given for one
// frame and reused by every frame are not (offset 0). setFrameViewsType6 (called once per reconstruction) works out the offsets, geomViewStartType6 and
// normViewStartType6 give them. The corrections of the measurement domain (corrVector) are offset by the frame in their own way (see type6Weights)
inline void setFrameViewsType6(const scalarStruct& inputScalars, Weighting& w_vec, const int64_t* pituus) {
	const size_t Nt = static_cast<size_t>(inputScalars.Nt);
	w_vec.frameStart6.assign(Nt + 1U, 0);
	w_vec.geomAllFrames6 = false;
	w_vec.normAllFrames6 = false;
	if (pituus == nullptr || Nt <= 1U)
		return;
	const size_t subsets = static_cast<size_t>(inputScalars.subsets);
	for (size_t tt = 0; tt <= Nt; tt++)
		w_vec.frameStart6[tt] = pituus[tt * subsets] - pituus[0];
	const int64_t total = w_vec.frameStart6[Nt];
	if (total <= w_vec.frameStart6[1])
		return;
	// The geometry holds all the frames when it has as many views as all the frames together (a geometry for one frame has the views of the frame)
	const int64_t views = w_vec.geomViews6 > 0 ? w_vec.geomViews6 : inputScalars.nProjections;
	w_vec.geomAllFrames6 = views == total;
	// The normalization has one image for each of the views of all the frames (a normalization with one image for each detector head is the same for every frame)
	const size_t nDet = static_cast<size_t>(inputScalars.nRowsD) * static_cast<size_t>(inputScalars.nColsD);
	w_vec.normAllFrames6 = inputScalars.normalization_correction && inputScalars.size_norm > 1ULL && nDet > 0ULL &&
		!(inputScalars.SPECT && inputScalars.normZ == inputScalars.nHeads) && inputScalars.size_norm == nDet * static_cast<size_t>(total);
}

// The offset of frame timestep in the geometry arrays and the DetectorVector (the first view of the frame there)
inline int64_t geomViewStartType6(const Weighting& w_vec, const int64_t timestep) {
	if (!w_vec.geomAllFrames6 || timestep <= 0 || static_cast<size_t>(timestep) + 1U >= w_vec.frameStart6.size())
		return 0;
	return w_vec.frameStart6[static_cast<size_t>(timestep)];
}

// The offset of frame timestep in the per-projection normalization
inline int64_t normViewStartType6(const Weighting& w_vec, const int64_t timestep) {
	if (!w_vec.normAllFrames6 || timestep <= 0 || static_cast<size_t>(timestep) + 1U >= w_vec.frameStart6.size())
		return 0;
	return w_vec.frameStart6[static_cast<size_t>(timestep)];
}

// Whether the geometry of the frames differs, for a geometry that holds all the frames: then the sensitivity image of the projector is not the same for every frame
// (a frame with other views than the first one, or the same views in another order, is another system)
inline bool geometryDiffersBetweenFramesType6(const Weighting& w_vec) {
	if (!w_vec.geomAllFrames6 || w_vec.frameStart6.size() < 3U)
		return false;
	const size_t Nt = w_vec.frameStart6.size() - 1U;
	const int64_t len0 = w_vec.frameStart6[1] - w_vec.frameStart6[0];
	for (size_t tt = 1U; tt < Nt; tt++) {
		const int64_t start = w_vec.frameStart6[tt];
		if (w_vec.frameStart6[tt + 1U] - start != len0)
			return true;
		const size_t n = static_cast<size_t>(len0);
		if (w_vec.swivelAngles != nullptr && std::memcmp(w_vec.swivelAngles, w_vec.swivelAngles + start, n * sizeof(float)) != 0)
			return true;
		if (w_vec.distFloat != nullptr && std::memcmp(w_vec.distFloat, w_vec.distFloat + start, n * sizeof(float)) != 0)
			return true;
		if (w_vec.distFloat2 != nullptr && std::memcmp(w_vec.distFloat2, w_vec.distFloat2 + start, n * sizeof(float)) != 0)
			return true;
		if (w_vec.distInt != nullptr && std::memcmp(w_vec.distInt, w_vec.distInt + start, n * sizeof(int32_t)) != 0)
			return true;
		if (w_vec.distInt2 != nullptr && std::memcmp(w_vec.distInt2, w_vec.distInt2 + start, n * sizeof(int32_t)) != 0)
			return true;
	}
	return false;
}

// The attenuation map of frame timestep of the host attenuation array of the type-6 projector. A dynamic scan can have a map for every frame (the vaimennus
// cell of Nt maps is concatenated by the host code into one array), each as large as the map of a static scan: the image of the main volume or, with
// multi-resolution volumes, the whole grid (the union bounding box of the volumes, see buildType6AttenuationVolume). A shorter array (a single map) is
// the map of every frame. Returns the pointer to the map of the frame and the index of the map (0 for a single map), by which the attenuation sources
// of the maps are kept apart in the cache of the projector (see AttCacheType6)
typedef struct AttFrameType6_ {
	const float* atten = nullptr;
	int frame = 0;
} AttFrameType6;

inline AttFrameType6 attenuationFrameType6(const scalarStruct& inputScalars, const float* atten, const int timestep) {
	AttFrameType6 out;
	out.atten = atten;
	if (atten == nullptr || timestep <= 0 || inputScalars.Nt <= 1U)
		return out;
	size_t frameElems = static_cast<size_t>(inputScalars.Nx[0]) * static_cast<size_t>(inputScalars.Ny[0]) * static_cast<size_t>(inputScalars.Nz[0]);
	if (inputScalars.nMultiVolumes > 0) {
		float bxMap = inputScalars.bx[0], byMap = inputScalars.by[0], bzMap = inputScalars.bz[0];
		float bxHi = inputScalars.bx[0] + inputScalars.Nx[0] * inputScalars.dx[0];
		float byHi = inputScalars.by[0] + inputScalars.Ny[0] * inputScalars.dy[0];
		float bzHi = inputScalars.bz[0] + inputScalars.Nz[0] * inputScalars.dz[0];
		for (int jj = 1; jj <= static_cast<int>(inputScalars.nMultiVolumes); jj++) {
			bxMap = (std::min)(bxMap, inputScalars.bx[jj]);
			byMap = (std::min)(byMap, inputScalars.by[jj]);
			bzMap = (std::min)(bzMap, inputScalars.bz[jj]);
			bxHi = (std::max)(bxHi, inputScalars.bx[jj] + inputScalars.Nx[jj] * inputScalars.dx[jj]);
			byHi = (std::max)(byHi, inputScalars.by[jj] + inputScalars.Ny[jj] * inputScalars.dy[jj]);
			bzHi = (std::max)(bzHi, inputScalars.bz[jj] + inputScalars.Nz[jj] * inputScalars.dz[jj]);
		}
		frameElems = static_cast<size_t>(std::lround((bxHi - bxMap) / inputScalars.dx[0])) * static_cast<size_t>(std::lround((byHi - byMap) / inputScalars.dy[0])) *
			static_cast<size_t>(std::lround((bzHi - bzMap) / inputScalars.dz[0]));
	}
	if (frameElems == 0ULL || inputScalars.size_atten < frameElems * static_cast<size_t>(inputScalars.Nt))
		return out;
	out.frame = (std::min)(timestep, static_cast<int>(inputScalars.Nt) - 1);
	out.atten = atten + frameElems * static_cast<size_t>(out.frame);
	return out;
}

#if defined(METAL)
// NOTE: multi-resolution volumes (nMultiVolumes > 0) are not yet supported/validated
// Builds multi-resolution volume ii's own attenuation image from the host-provided attenuation map.
// The map is expected to cover the union bounding box of every multi-resolution volume (0..
// nMultiVolumes), at volume 0's own (finest) voxel size -- i.e. exactly the same "single dense array
// covering the full, pre-split FOV" convention already used for options.x0 / options.referenceImage
// with multi-resolution volumes. With nMultiVolumes == 0 the map IS volume 0's own grid and this is
// a plain, unchanged reinterpretation of the host buffer (bit-identical to the pre-existing
// behavior). CTAttenuation must be true: type 6's rotate-and-integrate attenuation model needs a full
// image-domain map, unlike the precomputed LOR/sinogram-domain attenuation factors used elsewhere
// when CTAttenuation is false, which this projector cannot consume.
inline af::array buildType6AttenuationVolume(const scalarStruct& inputScalars, const int ii, const float* atten) {
	if (!(inputScalars.attenuation_correction && inputScalars.CTAttenuation && (atten != nullptr)))
		return af::array();
	if (inputScalars.nMultiVolumes == 0)
		return af::array(inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], atten);

	float bxMap = inputScalars.bx[0], byMap = inputScalars.by[0], bzMap = inputScalars.bz[0];
	float bxHi = inputScalars.bx[0] + inputScalars.Nx[0] * inputScalars.dx[0];
	float byHi = inputScalars.by[0] + inputScalars.Ny[0] * inputScalars.dy[0];
	float bzHi = inputScalars.bz[0] + inputScalars.Nz[0] * inputScalars.dz[0];
	for (int jj = 1; jj <= inputScalars.nMultiVolumes; jj++) {
		bxMap = (std::min)(bxMap, inputScalars.bx[jj]);
		byMap = (std::min)(byMap, inputScalars.by[jj]);
		bzMap = (std::min)(bzMap, inputScalars.bz[jj]);
		bxHi = (std::max)(bxHi, inputScalars.bx[jj] + inputScalars.Nx[jj] * inputScalars.dx[jj]);
		byHi = (std::max)(byHi, inputScalars.by[jj] + inputScalars.Ny[jj] * inputScalars.dy[jj]);
		bzHi = (std::max)(bzHi, inputScalars.bz[jj] + inputScalars.Nz[jj] * inputScalars.dz[jj]);
	}
	const int NxMap = static_cast<int>(std::lround((bxHi - bxMap) / inputScalars.dx[0]));
	const int NyMap = static_cast<int>(std::lround((byHi - byMap) / inputScalars.dy[0]));
	const int NzMap = static_cast<int>(std::lround((bzHi - bzMap) / inputScalars.dz[0]));
	const int nx = inputScalars.Nx[ii], ny = inputScalars.Ny[ii], nz = inputScalars.Nz[ii];
	// Identity: no multi-resolution split, the map already IS volume 0's own grid
	if (ii == 0 && inputScalars.Nx[0] == static_cast<uint32_t>(NxMap) && inputScalars.Ny[0] == static_cast<uint32_t>(NyMap) &&
		inputScalars.Nz[0] == static_cast<uint32_t>(NzMap))
		return af::array(NxMap, NyMap, NzMap, atten);
	// Summed-area table and box averages in double on the host (see resampleAttenuationType6Host), the float result uploaded once
	const std::vector<float> res = resampleAttenuationType6Host(atten, bxMap, byMap, bzMap, inputScalars.dx[0], inputScalars.dy[0], inputScalars.dz[0],
		NxMap, NyMap, NzMap, inputScalars.bx[ii], inputScalars.by[ii], inputScalars.bz[ii], inputScalars.dx[ii], inputScalars.dy[ii], inputScalars.dz[ii],
		nx, ny, nz, true);
	return af::array(nx, ny, nz, res.data());
}

// SPECT forward projection (projector type 6)
inline void forwardProjectionType6(af::array& fProj, const Weighting& w_vec, AF_im_vectors& vec, const scalarStruct& inputScalars,
	const int64_t length, const int64_t uu, ProjectorClass& proj, const int ii = 0, const float* atten = nullptr, const int timestep = 0) {
	if (!mainCenteredType6(inputScalars)) {
		mexPrint("An off-centre FOV is not supported with the SPECT projector type 6 on the Metal backend");
		return;
	}
	if (DEBUG || inputScalars.verbose >= 3) {
		af::sync();
		proj.tStartLocal = std::chrono::steady_clock::now();
	}
	if (DEBUG || inputScalars.verbose >= 3)
		mexPrint("Starting SPECT forward projection");
	int64_t u1 = uu + geomViewStartType6(w_vec, timestep);
	const af::array apuArr = ownedType6(af::moddims(vec.im_os[timestep][ii], inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii]));

	// Pre-transfer/resample the attenuation volume if it doesn't yet exist (once per call, not
	// per projection angle below)
	af::array attenBase = buildType6AttenuationVolume(inputScalars, ii, attenuationFrameType6(inputScalars, atten, timestep).atten);
	const bool useAtten = (inputScalars.attenuation_correction && inputScalars.CTAttenuation && (atten != nullptr));

	for (int kk = 0; kk < length; kk++) {
		af::array kuvaRot;

        // 1. Rotate the image
        kuvaRot = rotateHelperType6(apuArr, inputScalars, proj, -w_vec.swivelAngles[u1], ii);
        if (DEBUG || inputScalars.verbose > 2)
            mexPrint("Projector 6 FP step 1 complete");

        // 2. Translate the image
        if (w_vec.distInt2[u1] != 0)
            kuvaRot = translateHelperType6(kuvaRot, 0, -w_vec.distInt2[u1], 0);
        if (DEBUG || inputScalars.verbose > 2)
            mexPrint("Projector 6 FP step 2 complete");

        // 3. Process and apply attenuation image
        if (useAtten) {
            // 3.1. rotate attenuation map
            af::array attenuationImage = rotateHelperType6(attenBase, inputScalars, proj, -w_vec.swivelAngles[u1], ii);

            // 3.2 translate attenuation map, with the same shift as the emission volume above so that the
            // two stay aligned
            if (w_vec.distInt2[u1] != 0)
                attenuationImage = translateHelperType6(attenuationImage, 0, -w_vec.distInt2[u1], 0);

            // 3.3 accumulate attenuation
            attenuationImage = af::accum(attenuationImage, 0);

            // 3.4 multiply each accumulated element by voxel size and exponentiate attenuation
			attenuationImage = af::exp(-inputScalars.dx[ii] * attenuationImage);

            // 3.5 normalize attenuation map
            //af::array attenuationImageSum = af::sum(attenuationImage, 0) / inputScalars.Nx[0];
            //af::array ones = af::constant(1.0, attenuationImageSum.dims(), attenuationImageSum.type());
            //attenuationImageSum = af::select(attenuationImageSum == 0, ones, attenuationImageSum);
            //attenuationImage /= attenuationImageSum;

            // 3.6 pointwise multiply with image volume
            kuvaRot *= attenuationImage;
            af::eval(kuvaRot);
        }
        if (DEBUG || inputScalars.verbose > 2)
            mexPrint("Projector 6 FP step 3 complete");

        // 4. Convolve with detector PSF
        /* w_vec.distInt[u1] interpretation
            * A positive value indicates that the detector is inside the FOV
            * * After shift, set first w_vec.distInt[u1] elements through z-axis to zero
            * A negative value indicates that the detector is outside the FOV
            * * After shift, set last w_vec.distInt[u1] elements through z-axis to zero
        */
        af::array PSF = translateHelperType6(w_vec.gFilter, 0, 0, w_vec.distInt[u1]);
        PSF = PSF(af::span, af::span, af::seq(0, inputScalars.Nx[ii]-1));
        kuvaRot = af::reorder(kuvaRot, 2, 1, 0);
        kuvaRot = af::convolve2(kuvaRot, PSF);
        kuvaRot = af::reorder(kuvaRot, 2, 1, 0);
        //af::eval(kuvaRot);
        if (DEBUG || inputScalars.verbose > 2)
            mexPrint("Projector 6 FP step 4 complete");

        // 5. Sum
		//kuvaRot = af::reorder(kuvaRot, 2, 1, 0);
        kuvaRot = af::sum(kuvaRot, 0);
		kuvaRot = af::reorder(kuvaRot, 1, 2, 0);
        kuvaRot /= inputScalars.Nx[ii];
        if (DEBUG || inputScalars.verbose > 2)
            mexPrint("Projector 6 FP step 5 complete");

		fProj(af::span, af::span, kk) += kuvaRot.copy();
		u1++;
	}
	if (DEBUG || inputScalars.verbose >= 3) {
		proj.tEndLocal = std::chrono::steady_clock::now();
		const std::chrono::duration<double> tDiff = proj.tEndLocal - proj.tStartLocal;
		mexPrintBase("SPECT forward projection completed in %f seconds\n", tDiff);
		mexEval();
	}
}

// Backprojects a single projection angle. The PSF slice and the attenuation map depend only on the angle,
// so the caller builds them once and shares them between the measurement and the sensitivity image
inline af::array backProjectionType6Angle(const af::array& proj2D, const af::array& PSF, const af::array& attenuationImage,
    const bool useAtten, const Weighting& w_vec, const scalarStruct& inputScalars, ProjectorClass& proj, const int64_t u1, const int ii) {
    af::array kuvaRot = af::reorder(proj2D, 1, 0, 2); // Transpose the projection image

    // 1. Smear the input FP across the image volume
    kuvaRot = af::tile(kuvaRot, 1, 1, inputScalars.Nx[ii]); // Repeat through z-axis
    kuvaRot /= inputScalars.Nx[ii];

    // 2. Convolve with detector PSF. The forward projection applies the attenuation before the
    // convolution, so the adjoint has to apply the convolution first; the two do not commute
    kuvaRot = af::convolve2(kuvaRot, PSF);
    // Into the image volume layout, once, instead of transposing around each step below
    kuvaRot = reorder(kuvaRot, 2, 1, 0);

    // 3. Attenuation correction. Not evaluated here, the translation below forces evaluation on its own
    // through the indexed assignment, so an eval would only break the JIT fusion with the reorder above
    if (useAtten)
        kuvaRot *= attenuationImage;

    // 4. Translate the image
    if (w_vec.distInt2[u1] != 0)
        kuvaRot = translateHelperType6(kuvaRot, 0, w_vec.distInt2[u1], 0);
    af::eval(kuvaRot);

    // 5. Rotate the image
    kuvaRot = rotateHelperType6(kuvaRot, inputScalars, proj, w_vec.swivelAngles[u1], ii);
    af::eval(kuvaRot);
    return kuvaRot;
}

// Computes the rotation-based projector backprojection
// Currently this computes both the backprojection itself as well as the (optional) sensitivity image at the same time
inline af::array backProjectionType6Helper(const af::array &fProj, const Weighting& w_vec, const scalarStruct &inputScalars, ProjectorClass& proj, const int64_t length, uint32_t timestep, int64_t u1, const int ii, const float* atten = nullptr,
    const af::array* fProjSens = nullptr, af::array* outputSens = nullptr) {
    const af::array fp = ownedType6(af::moddims(fProj, inputScalars.nRowsD, inputScalars.nColsD, length));
    const bool compSens = (fProjSens != nullptr && outputSens != nullptr);
    af::array fpSens;
    if (compSens)
        fpSens = ownedType6(af::moddims(*fProjSens, inputScalars.nRowsD, inputScalars.nColsD, length));
	af::array output = af::constant(0.f, inputScalars.Nx[ii] * inputScalars.Ny[ii] * inputScalars.Nz[ii], 1);
    if (compSens)
        *outputSens = af::constant(0.f, inputScalars.Nx[ii] * inputScalars.Ny[ii] * inputScalars.Nz[ii], 1);

    const bool useAtten = (inputScalars.attenuation_correction && inputScalars.CTAttenuation && (atten != nullptr));
	// Pre-transfer/resample the attenuation image (once per call, not per projection angle below)
	af::array attenBase;
	if (useAtten)
		attenBase = buildType6AttenuationVolume(inputScalars, ii, attenuationFrameType6(inputScalars, atten, static_cast<int>(timestep)).atten);

    for (int kk = 0; kk < length; kk++) {
        /* N:=w_vec.distInt[u1] interpretation
            * A positive value indicates that the detector is inside the FOV
            * * After shift, set first N elements through z-axis to zero
            * A negative value indicates that the detector is outside the FOV
            * * After shift, set last N elements through z-axis to zero
        */
        af::array PSF = translateHelperType6(w_vec.gFilter, 0, 0, w_vec.distInt[u1]);
        PSF = PSF(af::span, af::span, af::seq(0, inputScalars.Nx[ii]-1));

        af::array attenuationImage;
        if (useAtten) {
            // Rotate the attenuation map, translate it with the same shift the forward projection applies,
            // accumulate along the projection direction and exponentiate
            attenuationImage = rotateHelperType6(attenBase, inputScalars, proj, -w_vec.swivelAngles[u1], ii);
            if (w_vec.distInt2[u1] != 0)
                attenuationImage = translateHelperType6(attenuationImage, 0, -w_vec.distInt2[u1], 0);
            attenuationImage = af::accum(attenuationImage, 0);
            attenuationImage = af::exp(-inputScalars.dx[ii] * attenuationImage);
        }

        const af::array slice = fp(af::span, af::span, kk);
        output += af::flat(backProjectionType6Angle(slice, PSF, attenuationImage, useAtten, w_vec, inputScalars, proj, u1, ii));
        if (compSens) {
            const af::array sliceSens = fpSens(af::span, af::span, kk);
            *outputSens += af::flat(backProjectionType6Angle(sliceSens, PSF, attenuationImage, useAtten, w_vec, inputScalars, proj, u1, ii));
        }
		u1++;
    }

    if (w_vec.maskBP != nullptr) { // Apply BP mask
        if (!applyMaskBP(output, compSens ? outputSens : nullptr, inputScalars, w_vec, ii, static_cast<int>(timestep)))
            return af::array();
    }

    return output;
}

// SPECT backprojection (projector type 6)
// sensIdx is the vec.Summ timestep slot to use (0 for a shared sensitivity image, or the real timestep
// when inputScalars.sensPerTimestep is set -- see reconstructionAF.h). allowSensCompute additionally
// gates the sensitivity computation itself: with a shared sensitivity image and Nt > 1, only the first
// timestep processed for a given subset (tt == t0) should (re)compute/finalize it; the caller passes
// false for every later timestep so it reuses the already-finalized image untouched instead of
// overwriting it with its own (otherwise identical) computation.
inline void backprojectionType6(af::array& fProj, const Weighting& w_vec, AF_im_vectors& vec,
	const scalarStruct& inputScalars, const int64_t length, const int64_t uu, ProjectorClass& proj, uint32_t timestep, const uint32_t osa_iter = 0, const uint32_t iter = 0,
	const uint8_t compute_norm_matrix = 0, const uint32_t iter0 = 0, const int ii = 0, const float* atten = nullptr, const uint32_t sensIdx = 0, const bool allowSensCompute = true) {
	if (!mainCenteredType6(inputScalars)) {
		mexPrint("An off-centre FOV is not supported with the SPECT projector type 6 on the Metal backend");
		return;
	}
	if (DEBUG || inputScalars.verbose >= 3)
		mexPrint("Starting SPECT backprojection");
	const bool compSens = allowSensCompute && ((iter == iter0 && compute_norm_matrix == 2) || compute_norm_matrix == 1);
	af::array sensProj, sensOut;
	if (compSens) {
		if (DEBUG || inputScalars.verbose >= 3)
			mexPrint("Computing sensitivity image");
		sensProj = af::constant(1.f, inputScalars.nRowsD, inputScalars.nColsD, length);
	}

    vec.rhs_os[timestep][ii] = backProjectionType6Helper(fProj, w_vec, inputScalars, proj, length, timestep, uu + geomViewStartType6(w_vec, static_cast<int64_t>(timestep)), ii, atten,
		compSens ? &sensProj : nullptr, compSens ? &sensOut : nullptr);
	vec.rhs_os[timestep][ii](vec.rhs_os[timestep][ii] < inputScalars.epps && vec.rhs_os[timestep][ii] >= 0.f) = inputScalars.epps;

	if (compSens) {
		if (compute_norm_matrix == 2) {
			vec.Summ[sensIdx][ii][osa_iter] = sensOut;
			vec.Summ[sensIdx][ii][osa_iter](vec.Summ[sensIdx][ii][osa_iter] < inputScalars.epps) = inputScalars.epps;
		} else {
			vec.Summ[sensIdx][ii][0] = sensOut;
			vec.Summ[sensIdx][ii][0](vec.Summ[sensIdx][ii][0] < inputScalars.epps) = inputScalars.epps;
		}
		if (DEBUG || inputScalars.verbose >= 3)
			mexPrint("Sensitivity image computed");
	}
	if (w_vec.maskBP != nullptr) {
		if (!resetMaskedBP(vec.rhs_os[timestep][ii], compSens ? &vec.Summ[sensIdx][ii][compute_norm_matrix == 2 ? osa_iter : 0] : nullptr, inputScalars, w_vec, ii, static_cast<int>(timestep)))
			return;
	}
	if (DEBUG || inputScalars.verbose >= 3)
		mexPrint("SPECT backprojection complete");
}
#else
// The depth-dependent collimator-detector response for a single projection, one filter per image plane
// The Gaussian CDRF is separable, so it is stored as an axial (fSize x 1 x Nx) and a transaxial
// (1 x fSize x Nx) filter bank. Convolving with the two of them in turn is both considerably faster
// than a single 2D convolution and needs a fraction of the memory. The precomputed bank, which may be
// an arbitrary user-defined filter, is kept as a single non-separable bank
// firstPlane is the first image plane that is in front of the collimator; the planes behind it are not
// visible to the detector and can be skipped entirely
// fused means that the response is computed inside the fused per-view kernels, see type6FusedAllowed, in
// which case only firstPlane and radius are set and no filter bank is built
// The kernels make the weights of the response of their planes themselves, except when the radius is larger
// than what they can handle (ProjectorClass::cdrfRMax) or the region is large (type6GlobalWeights): then the weights of one view and volume (a plane of taps
// for each image plane, see the cdrfWeights kernel) are in wAx and wTr, and the radius of each plane in rAx and
// rTr, once prepareFusedType6 has made them. They are made for the planes of the depth window (DepthWindowType6)
typedef struct PSFType6_ {
    af::array full, axial, transaxial;
    af::array wAx, wTr;
    // The radius of every plane for the two directions, int32 (see the cdrfWeights kernel)
    af::array rAx, rTr;
    bool separable = false;
    bool fused = false;
    int firstPlane = 0;
    int radius = 0;
} PSFType6;

// The attenuation that the axial kernel applies to a region of the grid (N0 x N1 x N2), from an array that may be larger than the
// region and need not be the factors: element (x, y, z) of the region is element (x + x0, y + y0, z + z0) of arr. With dx == 0 the
// array holds the attenuation factors, with dx > 0 the accumulated attenuation (the sum along the depth, without the voxel size, which
// is dx) and the factor is exp(-dx * value), computed by the kernel. The strides of the array are its dimensions
typedef struct AttRefType6_ {
	const af::array* arr = nullptr;
	int x0 = 0, y0 = 0, z0 = 0;
	float dx = 0.f;
	// With several views per launch (see type6BatchSize) the array is a stack of the attenuations of the views, one slice for each, and
	// this view's is the slice'th one, zStride planes along z (the third dimension of the slices) apart: the offset of the region
	// along z is z0 + slice * zStride. Zero for the array of a single view
	int slice = 0, zStride = 0;
} AttRefType6;

// The response is computed on-the-fly instead of taken from the precomputed bank. The precomputed bank is
// sampled at the voxel size of the main volume and can only be shifted by an integer number of voxels, so it
// is not applicable to the multi-resolution volumes, nor to a grid that is not the volume itself, that is an
// off-centre main volume, which is rotated into an axis-centered grid as well. Those always use the
// on-the-fly response, which the setup makes sure is available
inline bool onTheFlyPSFType6(const scalarStruct& inputScalars, const int ii) {
    return inputScalars.CDRFOnTheFly || ii > 0 || inputScalars.nMultiVolumes > 0 || !mainCenteredXYType6(inputScalars);
}

// Whether the fused per-view pipeline (kernels cdrfWeights, cdrfAxial, cdrfTransFP, cdrfTransBP
// and the detector resampling ones, a handful of launches per view instead of ~30 ArrayFire operations, with
// the results added directly into the projections and the image) is used for the on-the-fly
// response. It is available on the GPU backends only, the CPU backend and the precomputed (non-separable)
// filter bank keep the ArrayFire pipeline
#ifndef CPU
inline bool type6FusedAllowed(const scalarStruct& inputScalars, const int ii) {
    return onTheFlyPSFType6(inputScalars, ii);
}

// The fused kernels either make the weights of the response themselves in shared memory (cheapest when a view and
// volume has few voxels, where the launch overhead dominates and the weights are made once per group of planes) or
// read them from global memory, after cdrfWeights has made them once per view and volume (cheaper when the region
// is large, since every one of the many groups would otherwise make them again). A view and volume uses the
// weights in global memory when its number of voxels, the planes of the depth window times the rows and the axial rows
// of the region that is processed, is at least this threshold (see type6GlobalWeights): 1000000 voxels (a centred
// 256^3 volume is above it, the small multi-resolution volumes are below). The two give the same result up to the order
// of the sums
constexpr int64_t type6WeightsThreshold = 1000000;

// The attenuation of one batch of views (the views k0 ... k0 + nb - 1 of the call that starts at projection uu and has length views): the
// accumulated attenuation of every group (acc) and where the volumes read it (refs[view][volume], whose arr point to acc). sig is what the
// arrays depend on (see attBatchSigType6), which is compared before they are reused
typedef struct AttBatchType6_ {
	int64_t uu = 0, length = 0, k0 = 0;
	int nb = 0;
	std::vector<double> sig;
	std::vector<af::array> acc;
	std::vector<std::vector<AttRefType6>> refs;
	int64_t bytes = 0;
} AttBatchType6;

// The cache of a ProjectorClass (ProjectorClass::type6AttCache), so that it is released with the projector, which exists for one reconstruction
// only: nothing is left over from a reconstruction to the next one. The sources are rebuilt (and the batches dropped) when the key of the
// sources (the host pointer of the map, the dimensions and the lattices) changes. A dynamic scan with a map for every frame (see
// attenuationFrameType6) has a set of sources for every map, each built the first time that its frame is projected and kept, as the frames
// alternate within every subset: Nt sets, the memory of Nt sources, but no resampling after the first one of each frame
typedef struct AttCacheType6_ {
	// The attenuation sources of the attenuation groups (the projector that works on every volume of a view together), a set for every map,
	// their keys, and whether the source (of the map that upFrame is) is in its rotation image (rotSlotAttenuationType6 + group), which the
	// image keeps from a call to the next one. frame is the map of the latest call
	std::vector<std::vector<double>> srcKey;
	std::vector<std::vector<af::array>> src;
	std::vector<char> uploaded;
	int upFrame = 0;
	int frame = 0;
	// The attenuation sources of the volume-by-volume projector (a volume each, the volume ii of the map f is volSrc[f * (nMultiVolumes + 1) + ii]), and their keys
	std::vector<std::vector<double>> volKey;
	std::vector<af::array> volSrc;
	// The attenuation of the batches of views that the forward projection of the latest subset made, and the bytes of the arrays
	std::vector<std::shared_ptr<AttBatchType6>> batches;
	int64_t batchBytes = 0;
} AttCacheType6;

inline AttCacheType6& attCacheType6(ProjectorClass& proj) {
	if (!proj.type6AttCache)
		proj.type6AttCache = std::make_shared<AttCacheType6>();
	return *std::static_pointer_cast<AttCacheType6>(proj.type6AttCache);
}

// What a source of the attenuation depends on: the host map (the pointer, the lattice of the whole map) and the axis
inline std::vector<double> attMapKeyType6(const float* atten, const scalarStruct& inputScalars) {
	std::vector<double> k;
	k.push_back(static_cast<double>(reinterpret_cast<uintptr_t>(atten)));
	k.push_back(static_cast<double>(inputScalars.nMultiVolumes));
	k.push_back(mainCenteredXYType6(inputScalars) ? 1.0 : 0.0);
	float cx = 0.f, cy = 0.f, cz = 0.f;
	axisType6(inputScalars, cx, cy, cz);
	k.push_back(cx); k.push_back(cy); k.push_back(cz);
	for (int ii = 0; ii <= static_cast<int>(inputScalars.nMultiVolumes); ii++) {
		k.push_back(static_cast<double>(inputScalars.Nx[ii])); k.push_back(static_cast<double>(inputScalars.Ny[ii])); k.push_back(static_cast<double>(inputScalars.Nz[ii]));
		k.push_back(inputScalars.bx[ii]); k.push_back(inputScalars.by[ii]); k.push_back(inputScalars.bz[ii]);
		k.push_back(inputScalars.dx[ii]); k.push_back(inputScalars.dy[ii]); k.push_back(inputScalars.dz[ii]);
	}
	return k;
}

// The key of the attenuation sources of the groups: the map and the lattice of every group
inline std::vector<double> attGroupKeyType6(const float* atten, const scalarStruct& inputScalars, const std::vector<AttGroupType6>& groups) {
	std::vector<double> k = attMapKeyType6(atten, inputScalars);
	k.push_back(static_cast<double>(groups.size()));
	for (size_t gg = 0; gg < groups.size(); gg++) {
		const AttGroupType6& grp = groups[gg];
		k.push_back(static_cast<double>(grp.members.size()));
		for (size_t mm = 0; mm < grp.members.size(); mm++)
			k.push_back(grp.members[mm]);
		k.push_back(grp.attNx); k.push_back(grp.srcNx); k.push_back(grp.srcNy); k.push_back(grp.rows);
		k.push_back(grp.bz); k.push_back(grp.nz); k.push_back(grp.padLo); k.push_back(grp.padHi);
		for (size_t zz = 0; zz < grp.zOff.size(); zz++)
			k.push_back(grp.zOff[zz]);
	}
	return k;
}

// The key of the attenuation source of volume ii of the volume-by-volume projector: the map and the lattice of the volume
inline std::vector<double> attVolumeKeyType6(const float* atten, const scalarStruct& inputScalars, const int ii) {
	std::vector<double> k = attMapKeyType6(atten, inputScalars);
	int gridNx = 0, gridNy = 0, attNx = 0, srcNx = 0, srcNy = 0;
	rotGridType6(inputScalars, ii, gridNx, gridNy);
	attGridType6(inputScalars, ii, gridNx, gridNy, attNx, srcNx, srcNy);
	k.push_back(ii); k.push_back(attNx); k.push_back(srcNx); k.push_back(srcNy);
	return k;
}

// The device memory that the retained attenuation of a subset may use: the fraction of the available memory that the option type6AttenuationCache
// gives (default 0.4, clamped to [0, 0.9]; 0 keeps nothing, so the backprojection recomputes the attenuation). Available memory is, for CUDA and HIP, the
// free memory of the device and the part of ArrayFire's pool that is not in use, and, for OpenCL, the whole memory of the device (which has no
// query for the free part)
inline int64_t attBudgetType6(ProjectorClass& proj, const scalarStruct& inputScalars) {
	double frac = static_cast<double>(inputScalars.type6AttenuationCache);
	if (!(frac > 0.0))
		return 0;
	frac = (std::min)(frac, 0.9);
	int64_t avail = proj.getGlobalMem();
#if defined(CUDA) || defined(HIP)
	size_t memFree = 0ULL, memTotal = 0ULL;
	if (cuMemGetInfo(&memFree, &memTotal) == 0) {
		size_t allocBytes = 0, allocBuffers = 0, lockBytes = 0, lockBuffers = 0;
		af::deviceMemInfo(&allocBytes, &allocBuffers, &lockBytes, &lockBuffers);
		avail = static_cast<int64_t>(memFree) + static_cast<int64_t>(allocBytes > lockBytes ? allocBytes - lockBytes : 0);
	}
#endif
	return (std::max)(static_cast<int64_t>(frac * static_cast<double>(avail)), static_cast<int64_t>(0));
}


// The number of views that the fused pipeline processes per kernel launch (forwardProjectionType6AllBatched and
// backprojectionType6AllBatched): the environment variable TYPE6_BATCH, at most 64. The default (variable not set or invalid) is
// 8. 1 gives the pipeline that processes one view at a time (the previous behaviour). B > 1 is bitwise identical to B = 1 (verified on
// CUDA) and only batches views on the in-kernel-weights path. It is read once, when the projector is first used
inline int type6BatchSize() {
    static const int batch = []() {
        int b = 8;
#ifdef _WIN32
        char* val = nullptr;
        size_t len = 0;
        if (_dupenv_s(&val, &len, "TYPE6_BATCH") == 0 && val != nullptr) {
            const long v = std::strtol(val, nullptr, 10);
            if (v >= 1)
                b = static_cast<int>((std::min)(v, 64L));
            free(val);
        }
#else
        const char* val = std::getenv("TYPE6_BATCH");
        if (val != nullptr) {
            const long v = std::strtol(val, nullptr, 10);
            if (v >= 1)
                b = static_cast<int>((std::min)(v, 64L));
        }
#endif
        return b;
    }();
    return batch;
}

// The largest number of elements that one of the stacked arrays of a batch of views may have (the kernels index them with 32 bits and
// a batch is not worth a very large array)
constexpr int64_t type6BatchMaxElements = static_cast<int64_t>(1) << 25;

// The number of views per launch that the projectors of the volumes use (type6BatchSize), 1 meaning that they process the views one by
// one. It is 1 unless every volume can use the fused pipeline, and it is limited so that the stacked arrays of a batch, of the largest
// rotated grid and of the largest attenuation group, stay below type6BatchMaxElements and the launches below 65535 groups along their
// third dimension (the batched kernels, see the end of the type-6 projectors below, put the views there)
inline int type6BatchSizeEffective(const scalarStruct& inputScalars) {
    int batch = type6BatchSize();
    if (batch <= 1)
        return 1;
    int64_t slice = 1;
    int zMax = 1;
    for (int ii = 0; ii <= static_cast<int>(inputScalars.nMultiVolumes); ii++) {
        if (!type6FusedAllowed(inputScalars, ii))
            return 1;
        int gNx = 0, gNy = 0;
        rotGridType6(inputScalars, ii, gNx, gNy);
        const int gNz = static_cast<int>(inputScalars.Nz[ii]) + padZLoType6(inputScalars, ii) + padZHiType6(inputScalars, ii);
        slice = (std::max)(slice, static_cast<int64_t>(roundUpType6(gNx, type6AllocQuantum)) * roundUpType6(gNy, type6AllocQuantum) * roundUpType6(gNz, type6AllocQuantum));
        zMax = (std::max)(zMax, gNz);
    }
    const std::vector<AttGroupType6> groups = attGroupsType6(inputScalars);
    for (size_t gg = 0; gg < groups.size(); gg++) {
        const int gNz = groups[gg].padLo + groups[gg].nz + groups[gg].padHi;
        slice = (std::max)(slice, static_cast<int64_t>(roundUpType6(groups[gg].attNx, type6AllocQuantum)) * roundUpType6(groups[gg].rows, type6AllocQuantum)
            * roundUpType6(gNz, type6AllocQuantum));
        zMax = (std::max)(zMax, gNz);
    }
    batch = static_cast<int>((std::min)(static_cast<int64_t>(batch), type6BatchMaxElements / slice));
    batch = (std::min)(batch, 65535 / zMax);
    return (std::max)(batch, 1);
}

#else
inline bool type6FusedAllowed(const scalarStruct&, const int) {
    return false;
}
#endif

// The grid that one volume is rotated into for one projection. Without multi-resolution the volume is
// rotated into itself, which is what the projector has always done. A multi-resolution subvolume is
// neither centered on the rotation axis of the gantry nor sampled at the pitch of the detector, so it
// is rotated into a grid that is centered on the axis instead and keeps the voxel size of the
// subvolume. Its projection then lands on a coarser lattice than the detector and is resampled onto it
typedef struct GridType6_ {
    // Dimensions of the grid. The axial direction is not rotated, so it is always that of the volume
    int Nx = 0, Ny = 0, Nz = 0;
    // Distance from the first plane of the grid to the collimator surface, in voxels of the grid
    float dist = 0.f;
    // Lateral shift of the detector panel, in voxels of the grid
    float panelShift = 0.f;
    // Position of the rotation axis with respect to the center of the volume, in voxels of the volume.
    // Zero for the main volume, which is centered on the axis
    float axisX = 0.f, axisY = 0.f;
    // Whether the projection of the grid is on the detector lattice or has to be resampled onto it
    bool resample = false;
    // The variance (voxels^2) that the rotation of the view adds transaxially (rotVarType6, plus that of the resampling onto the detector
    // with resample), which the transaxial response is narrowed by. Negative: the constant sigmaCoeff^2 of before (TYPE6_ROTCOMP=0, Metal)
    float rotVar = -1.f;
} GridType6;

inline void offsetsFPType6(const Weighting& w_vec, const GridType6& grid, const int64_t u1, float& offX, float& offY);

inline GridType6 gridType6(const Weighting& w_vec, const scalarStruct& inputScalars, const int64_t u1, const int ii = 0) {
    GridType6 grid;
    rotGridType6(inputScalars, ii, grid.Nx, grid.Ny);
    // With multi-resolution the grid is padded axially, see computePadZType6. Grid slice k is slice
    // k - padZLoType6 of the volume, and the grid is padded only at the ends that the detector extends past
    grid.Nz = static_cast<int>(inputScalars.Nz[ii]) + padZLoType6(inputScalars, ii) + padZHiType6(inputScalars, ii);
    // The main volume is on the lattice of the detector, but a detector plane that is wider than the grid (an
    // uncropped projection padded to the extended FOV) is handled by the same resampling. With ii == 0 that
    // is an exact placement of the grid onto the plane, the tent weights are 0 or 1 at integer offsets
    // The detector is centered on the rotation axis and on z = 0, so a main volume that is not centered on them
    // (an object offset) is not centered on the detector either and needs the exact mapping as well
    // With different padding at the two ends the main volume is not centered in its grid either, which the non-resampled placement assumes
    grid.resample = ii > 0 || !mainCenteredType6(inputScalars) || padZLoType6(inputScalars, ii) != padZHiType6(inputScalars, ii) || static_cast<int>(inputScalars.nRowsD) > grid.Ny || static_cast<int>(inputScalars.nColsD) > grid.Nz;
    // The geometry is given in voxels of the main volume, so it is converted into voxels of this grid.
    // The collimator surface is in the same place for every volume, but the grids start at different
    // distances from the rotation axis
    const float collimator = (w_vec.distFloat[u1] - .5f * static_cast<float>(inputScalars.Nx[0])) * inputScalars.dx[0];
    grid.dist = (collimator + .5f * static_cast<float>(grid.Nx) * inputScalars.dx[ii]) / inputScalars.dx[ii];
    grid.panelShift = w_vec.distFloat2[u1] * inputScalars.dx[0] / inputScalars.dx[ii];
    // The detector window of the placement on the lattice (cropType6) starts at row (Ny - nRowsD) / 2 (integer division), which is half a row off the rotation
    // axis, the exact centre of the grid (see the rotation kernels), when Ny - nRowsD is odd. The half row is made up by the lateral shift of the panel
    if (!grid.resample && ((grid.Ny - static_cast<int>(inputScalars.nRowsD)) % 2) != 0)
        grid.panelShift += 0.5f;
    if (ii > 0 || !mainCenteredXYType6(inputScalars)) {
        // The rotation axis is the origin. The main volume is rotated into itself when it is centered on it
        float cx = 0.f, cy = 0.f, cz = 0.f;
        axisType6(inputScalars, cx, cy, cz);
        grid.axisX = (cx - (inputScalars.bx[ii] + .5f * static_cast<float>(inputScalars.Nx[ii]) * inputScalars.dx[ii])) / inputScalars.dx[ii];
        grid.axisY = (cy - (inputScalars.by[ii] + .5f * static_cast<float>(inputScalars.Ny[ii]) * inputScalars.dy[ii])) / inputScalars.dy[ii];
    }
#ifndef METAL
    if (type6RotComp()) {
        // The transaxial variance of the interpolation of the rotation of this view (see rotVarType6), with the very values that the rotation
        // kernel is given (type6PutRotation, rotateCustomAF): the rotation by -swivel, the offsets of offsetsFPType6, the volume rotated
        // into the grid. The backprojection and the sensitivity image use the same grid, so the response of the two directions stays the same
        float offX = 0.f, offY = 0.f;
        offsetsFPType6(w_vec, grid, u1, offX, offY);
        const float angleRad = (-w_vec.swivelAngles[u1]) * static_cast<float>(M_PI) / 180.f;
        const float cosa = std::cos(-angleRad);
        const float sina = std::sin(-angleRad);
        double var = rotVarType6(static_cast<double>(cosa), static_cast<double>(sina), static_cast<double>(offX), static_cast<double>(offY),
            static_cast<int>(inputScalars.Nx[ii]), static_cast<int>(inputScalars.Ny[ii]), grid.Nx, grid.Ny);
        // With the projection resampled onto the detector (a map that is the same for every view up to a whole number of rows, see mapType6)
        // that interpolation blurs it as well
        if (grid.resample) {
            const double trA = static_cast<double>(inputScalars.dy[0]) / static_cast<double>(inputScalars.dy[ii]);
            const double trB = .5 * static_cast<double>(grid.Ny) + (.5 - .5 * static_cast<double>(inputScalars.nRowsD)) * trA - .5;
            var += detVarType6Compute(static_cast<int>(inputScalars.nRowsD), trA, trB);
        }
        grid.rotVar = static_cast<float>(var);
    }
#endif
    return grid;
}

// The projection image does not necessarily cover the whole FOV. With cropped projections it covers
// only the detector panel, which the rotation places at the center of the image. In that case just the
// part of the volume that can reach those detector bins has to be convolved, that is the detector
// window widened by the radius of the CDRF. The result is identical to convolving the whole volume and
// cropping afterwards, since the filters are truncated at that radius. The same applies to the
// multi-resolution volumes, whose grid mostly reaches past the detector at any given angle
typedef struct CropType6_ {
    // Region of the convolved volume that is computed, dimension 0 is axial (Nz) and 1 transaxial (Ny)
    int axialStart = 0, axialLength = 0, transStart = 0, transLength = 0;
    // Position of the detector window inside that region, used when the region is on the detector
    // lattice, i.e. without multi-resolution
    int axialOffset = 0, transOffset = 0;
    bool cropped = false;
    // Nothing in the volume can reach the detector, which happens with the multi-resolution volumes
    // that extend past the axial coverage of the detector
    bool empty = false;
} CropType6;

inline CropType6 cropType6(const scalarStruct& inputScalars, const PSFType6& psf, const GridType6& grid, const int ii = 0) {
    const int Ny = grid.Ny;
    const int Nz = grid.Nz;
    const int nRows = static_cast<int>(inputScalars.nRowsD); // Transaxial detector bins
    const int nCols = static_cast<int>(inputScalars.nColsD); // Axial detector bins
    CropType6 crop;
    // The axial margin of the window is the axial radius of the response, which the padding of the grid
    // (padZType6, the bound of the padding) is made for, so with multi-resolution it is at most that. The margin is the same at both
    // ends and is clamped to the grid, which is padded only as far as the detector extends past the volume, so the window stays
    // valid: where the grid has less padding the window just ends at the edge of the grid. The radius of the response is
    // shared by both directions and can be larger because of the transaxial one
    const int padZ = padZType6(inputScalars, ii);
    const int padLo = padZLoType6(inputScalars, ii);
    const int marginZ = usesPadZType6(inputScalars) ? (std::min)(psf.radius, padZ) : psf.radius;
    if (!grid.resample) {
        crop.axialLength = Nz;
        crop.transLength = Ny;
        // A projection larger than the FOV cannot be cropped. That is not a supported combination and the
        // dimension mismatch is reported by the assignment further down, but the window has to stay valid
        crop.cropped = (nRows < Ny || nCols < Nz) && nRows <= Ny && nCols <= Nz;
        if (!crop.cropped)
            return crop;
        // The detector window is centered on the rotation axis
        const int z0 = (Nz - nCols) / 2;
        const int y0 = (Ny - nRows) / 2;
        crop.axialStart = (std::max)(z0 - marginZ, 0);
        crop.transStart = (std::max)(y0 - psf.radius, 0);
        crop.axialLength = (std::min)(z0 + nCols + marginZ, Nz) - crop.axialStart;
        crop.transLength = (std::min)(y0 + nRows + psf.radius, Ny) - crop.transStart;
        crop.axialOffset = z0 - crop.axialStart;
        crop.transOffset = y0 - crop.transStart;
        return crop;
    }
    // With multi-resolution the detector panel is, in general, neither centered in the grid nor an
    // integer number of voxels wide, so the window is worked out from the positions of the outermost
    // detector bins. Transaxially the panel is centered on the rotation axis, which the rotation places
    // at the center of the grid, and axially on the center of the main volume. The window has to hold
    // the two voxels that the resampling interpolates between, hence the extra voxel at the far end
    const float halfY = (.5f * static_cast<float>(nRows) - .5f) * inputScalars.dy[0] / inputScalars.dy[ii];
    // The detector is centered axially on z = 0 (see axisType6)
    float axisXc = 0.f, axisYc = 0.f, centerZ = 0.f;
    axisType6(inputScalars, axisXc, axisYc, centerZ);
    const float halfZ = (.5f * static_cast<float>(nCols) - .5f) * inputScalars.dz[0] / inputScalars.dz[ii];
    const float panelY = .5f * static_cast<float>(Ny) - .5f;
    // The volume starts padLo slices into the (axially padded) grid
    const float panelZ = (centerZ - inputScalars.bz[ii]) / inputScalars.dz[ii] - .5f + static_cast<float>(padLo);
    const int y0 = static_cast<int>(std::floor(panelY - halfY));
    const int y1 = static_cast<int>(std::floor(panelY + halfY)) + 2;
    const int z0 = static_cast<int>(std::floor(panelZ - halfZ));
    const int z1 = static_cast<int>(std::floor(panelZ + halfZ)) + 2;
    crop.transStart = (std::max)(y0 - psf.radius, 0);
    crop.axialStart = (std::max)(z0 - marginZ, 0);
    crop.transLength = (std::min)(y1 + psf.radius, Ny) - crop.transStart;
    crop.axialLength = (std::min)(z1 + marginZ, Nz) - crop.axialStart;
    crop.cropped = true;
    crop.empty = crop.transLength <= 0 || crop.axialLength <= 0;
    return crop;
}

// Maps the detector onto the cropped region of the rotated grid. Detector column c and row r fall on
// the position axA * c + axB along the axial and trA * r + trB along the transaxial direction of the
// region. Without multi-resolution the two lattices are the same and the map is a plain shift
typedef struct MapType6_ {
    float axA = 1.f, axB = 0.f, trA = 1.f, trB = 0.f;
} MapType6;

inline MapType6 mapType6(const scalarStruct& inputScalars, const CropType6& crop, const GridType6& grid, const int ii = 0) {
    MapType6 map;
    if (!grid.resample) {
        map.axB = static_cast<float>(crop.axialOffset);
        map.trB = static_cast<float>(crop.transOffset);
        return map;
    }
    float axisXc = 0.f, axisYc = 0.f, centerZ = 0.f;
    axisType6(inputScalars, axisXc, axisYc, centerZ);
    map.axA = inputScalars.dz[0] / inputScalars.dz[ii];
    map.axB = (centerZ - inputScalars.bz[ii]) / inputScalars.dz[ii]
        + (.5f - .5f * static_cast<float>(inputScalars.nColsD)) * map.axA - .5f + static_cast<float>(padZLoType6(inputScalars, ii)) - static_cast<float>(crop.axialStart);
    map.trA = inputScalars.dy[0] / inputScalars.dy[ii];
    map.trB = .5f * static_cast<float>(grid.Ny)
        + (.5f - .5f * static_cast<float>(inputScalars.nRowsD)) * map.trA - .5f - static_cast<float>(crop.transStart);
    return map;
}

// The rotation writes only the cropped region, which is what the forward projection needs. With xEnd > xOff
// only the planes xOff ... xEnd - 1 of the grid are written (the depth window of the fused pipeline)
inline RotWindow rotWindowOut(const CropType6& crop, const GridType6& grid, const scalarStruct& inputScalars, const int ii = 0, const int xOff = 0, const int xEnd = 0) {
    RotWindow win;
    if (xEnd > xOff) {
        win.outX = xOff;
        win.outNxW = xEnd - xOff;
    }
    win.srcNx = static_cast<int>(inputScalars.Nx[ii]);
    win.srcNyFull = static_cast<int>(inputScalars.Ny[ii]);
    // The volume is placed padLo slices into the axially padded grid. The rotation reads the source at the
    // slice it writes, so the source dimension is that of the grid and the volume is the band of it
    const int padZ = padZLoType6(inputScalars, ii);
    win.srcNzFull = grid.Nz;
    win.outNx = grid.Nx; win.outNyFull = grid.Ny;
    win.outY = crop.transStart; win.outNy = crop.transLength;
    win.outZ = crop.axialStart; win.outNz = crop.axialLength;
    win.srcY = 0; win.srcNy = win.srcNyFull;
    win.srcZ = padZ; win.srcNz = static_cast<int>(inputScalars.Nz[ii]);
    return win;
}
// The rotation reads only the cropped region and writes the whole volume, which is what the
// backprojection needs to turn the cropped result back into an image. With xEnd > xOff the source holds
// only the planes xOff ... xEnd - 1 of the grid (the depth window of the fused pipeline)
inline RotWindow rotWindowIn(const CropType6& crop, const GridType6& grid, const scalarStruct& inputScalars, const int ii = 0, const int xOff = 0, const int xEnd = 0) {
    RotWindow win;
    if (xEnd > xOff) {
        win.srcX = xOff;
        win.srcNxB = xEnd - xOff;
    }
    win.srcNx = grid.Nx; win.srcNyFull = grid.Ny; win.srcNzFull = grid.Nz;
    win.outNx = static_cast<int>(inputScalars.Nx[ii]);
    win.outNyFull = static_cast<int>(inputScalars.Ny[ii]);
    win.outY = 0; win.outNy = win.outNyFull;
    // The volume is the part of the axially padded grid that starts padLo slices in. The rotation reads the
    // source at slice outZ + k to write slice k of the output
    win.outZ = padZLoType6(inputScalars, ii); win.outNz = static_cast<int>(inputScalars.Nz[ii]);
    win.srcY = crop.transStart; win.srcNy = crop.transLength;
    win.srcZ = crop.axialStart; win.srcNz = crop.axialLength;
    return win;
}
// The detector window inside the cropped region
inline af::seq windowAxialSeq(const CropType6& crop, const scalarStruct& inputScalars) { return af::seq(crop.axialOffset, crop.axialOffset + static_cast<int>(inputScalars.nColsD) - 1); }
inline af::seq windowTransSeq(const CropType6& crop, const scalarStruct& inputScalars) { return af::seq(crop.transOffset, crop.transOffset + static_cast<int>(inputScalars.nRowsD) - 1); }

// Linear (tent) interpolation weights from the nSrc samples of a lattice onto nDet detector positions,
// as a (nDet x nSrc) matrix. Detector position d samples the lattice at A * d + B. The samples that
// are missing outside the lattice have no weight, i.e. the projection reads as zero there
inline af::array detectorWeightsType6(const int nDet, const int nSrc, const float A, const float B) {
    const af::array d = af::range(af::dim4(nDet, nSrc), 0);
    const af::array s = af::range(af::dim4(nDet, nSrc), 1);
    return (af::max)(1.f - af::abs(d * A + B - s), 0.f);
}

// Resamples the projection of a multi-resolution volume, which is on the lattice of its rotated grid
// (axial x transaxial), onto the detector (nColsD x nRowsD). The projection is the line integral through
// the volume, so the value of a detector bin is the projection at the position of that bin. The
// interpolation is separable, so it is exactly Wax * P * Wtr^T with the weights above
inline af::array toDetectorType6(const af::array& proj, const CropType6& crop, const MapType6& map, const scalarStruct& inputScalars) {
    const af::array Wax = detectorWeightsType6(static_cast<int>(inputScalars.nColsD), crop.axialLength, map.axA, map.axB);
    const af::array Wtr = detectorWeightsType6(static_cast<int>(inputScalars.nRowsD), crop.transLength, map.trA, map.trB);
    return af::matmul(af::matmul(Wax, proj), Wtr, AF_MAT_NONE, AF_MAT_TRANS);
}

// The exact adjoint of the above, Wax^T * M * Wtr, from the measurements (nColsD x nRowsD) onto the
// lattice of the rotated grid (axial x transaxial). There is no area factor, the transpose of the
// sampling is what the adjoint of the forward projection is
inline af::array fromDetectorType6(const af::array& meas, const CropType6& crop, const MapType6& map, const scalarStruct& inputScalars) {
    const af::array Wax = detectorWeightsType6(static_cast<int>(inputScalars.nColsD), crop.axialLength, map.axA, map.axB);
    const af::array Wtr = detectorWeightsType6(static_cast<int>(inputScalars.nRowsD), crop.transLength, map.trA, map.trB);
    return af::matmul(af::matmul(Wax, meas, AF_MAT_TRANS, AF_MAT_NONE), Wtr);
}

// The (not yet normalized) response of the square model (CDRFModel = 'square') of the planes, for the ArrayFire pipeline (the CPU backend, the
// GPU backends make the weights in the kernels, see cdrfTap in auxKernels.cl, whose formulas these are, in double precision): the
// triangle of FWHM (half base) Rg, convolved with the intrinsic Gaussian of FWHM C, point sampled at the taps -r ... r. Rg (mm) is
// 1 x 1 x nPlanes, dd the voxel size of the direction (mm), trans whether the rotation blur is taken out of the response. The result is
// (2 r + 1) x 1 x nPlanes, float
inline af::array squareBankType6(const af::array& Rg, const float C, const float dd, const float sigmaCoeff, const bool trans, const int r, const int nPlanes) {
    const double sMin = 1e-3;
    const double si = static_cast<double>(C) / dd * sigmaCoeff;
    const double comp = trans ? static_cast<double>(sigmaCoeff) * sigmaCoeff : 0.;
    const double v = si * si - comp;
    const af::array Rv = Rg.as(f64) / static_cast<double>(dd);
    af::array R, s;
    if (v >= sMin * sMin) {
        R = Rv;
        s = af::constant(std::sqrt(v), Rv.dims(), f64);
    }
    else {
        // The intrinsic Gaussian cannot take out the variance of the rotation blur, so the rest is taken out of the triangle (variance R^2 / 6)
        R = af::sqrt((af::max)(Rv * Rv - 6. * (sMin * sMin - v), 0.));
        s = af::constant(sMin, Rv.dims(), f64);
    }
    const int fSize = 2 * r + 1;
    const af::array u = af::tile(af::abs(af::range(af::dim4(fSize), 0, f64) - static_cast<double>(r)), 1, 1, nPlanes);
    const af::array Rt = af::tile(R, fSize, 1, 1);
    const af::array St = af::tile(s, fSize, 1, 1);
    // P(x) = E[(S - x)+], S ~ N(0, s^2)
    auto P = [&](const af::array& x) {
        const af::array z = x / St;
        return St * af::exp(-0.5 * z * z) * 0.3989422804014327 - x * 0.5 * af::erfc(z * 0.7071067811865476);
    };
    const af::array closed = (P(u + Rt) - 2. * P(u) + P(af::abs(u - Rt)) + af::select(u < Rt, Rt - u, af::constant(0., u.dims(), f64))) / (Rt * Rt);
    // R small with respect to s: the series in the even moments of the triangle, as in cdrfSquareTap
    const af::array a = Rt / St;
    const af::array z2 = (u / St) * (u / St);
    const af::array a2 = a * a;
    const af::array he2 = z2 - 1.;
    const af::array he4 = (z2 - 6.) * z2 + 3.;
    const af::array he6 = ((z2 - 15.) * z2 + 45.) * z2 - 15.;
    const af::array he8 = (((z2 - 28.) * z2 + 210.) * z2 - 420.) * z2 + 105.;
    const af::array f = 1. + a2 * (he2 / 12. + a2 * (he4 / 360. + a2 * (he6 / 20160. + a2 * he8 / 1814400.)));
    const af::array series = af::exp(-0.5 * z2) * 0.3989422804014327 / St * f;
    return af::select(a < 0.5, series, closed).as(f32);
}

// The filters are computed on-the-fly, instead of using a precomputed bank, because with swivelling
// detector heads the distance from the collimator to the image planes is different for every projection
// and, with up to several thousand projections, precomputing a separate bank for each of them is not
// feasible. Shifting a single precomputed bank, on the other hand, is limited to integer voxel steps,
// which is accurate enough only when the detector distance is identical for every projection
// The resolution is FWHM_axis(z) = sqrt((A_axis * z + B_axis)^2 + C^2), z the distance from the
// collimator surface in mm, with a separate (A, B) for the axial and for the transaxial direction (see
// scalarStruct::CDRFCoeffs6). The axial FWHM is converted to voxels with dz, the transaxial one with dy
// and the depth of a plane with dx
// With bank set to false only firstPlane and radius are worked out (on the host), which is all that the
// fused kernels need since they compute the weights themselves
inline PSFType6 psfHelperType6(const Weighting& w_vec, const scalarStruct& inputScalars, const GridType6& grid, const int ii = 0, const bool bank = true) {
    const float dx = inputScalars.dx[ii];
    const float dy = inputScalars.dy[ii];
    const float dz = inputScalars.dz[ii];
    const float Az = inputScalars.CDRFCoeffs6[0], Bz = inputScalars.CDRFCoeffs6[1];
    const float Axy = inputScalars.CDRFCoeffs6[2], Bxy = inputScalars.CDRFCoeffs6[3];
    const float C = inputScalars.CDRFCoeffs6[4];
    // The response is applied to the rotated grid, so both the number of planes and the distance to the
    // collimator are those of the grid and not those of the volume
    const int Nx = grid.Nx;
    const float dist = grid.dist;
    const float sigmaCoeff = 1.f / (2.f * std::sqrt(2.f * std::log(2.f)));
    PSFType6 psf;
    psf.separable = true;

    // Planes behind the collimator surface are not visible to the detector. Since the image is rotated
    // so that the detector is parallel with it, these are always the first planes.
    // The start is rounded down to a multiple of planeStep and the filter radius up to a multiple of
    // radiusStep further down. Without that the array shapes would be different for almost every
    // projection, ArrayFire could never reuse a pooled buffer, and both the allocation traffic and the
    // size of the pool would grow badly. The planes that the rounding adds are masked out below, so the
    // result does not change
    const int planeStep = 16;
    psf.firstPlane = (std::min)((std::max)(static_cast<int>(std::ceil(dist - 0.5f)), 0), Nx - 1);
    psf.firstPlane = (psf.firstPlane / planeStep) * planeStep;
    const int nPlanes = Nx - psf.firstPlane;

    // Truncate the Gaussians where they fall below 1e-6 of their peak value. The widest one is at one
    // of the two ends of the plane range, since the resolution is monotonic in the distance, so the
    // radius is worked out on the host. Reading it back from the device instead would block the queue
    // once per projection and leave the device idle for the rest of the loop. The radius is shared by
    // the two directions, so it is that of the widest of them (without the compensation above, which
    // only narrows the transaxial one)
    const float zFirst = (0.5f + static_cast<float>(psf.firstPlane) - dist) * dx;
    const float zLast = (static_cast<float>(Nx) - 0.5f - dist) * dx;
    // (With the square model, CDRFModel, the response reaches the half base of the triangle plus 5.2565 sigma of the intrinsic Gaussian)
    float extentMax = 0.f;
    for (int ee = 0; ee < 2; ee++) {
        const float zz = ee == 0 ? zFirst : zLast;
        const float RgEndZ = (std::max)(Az * zz + Bz, 0.f);
        const float RgEndXY = (std::max)(Axy * zz + Bxy, 0.f);
        const float eEndZ = cdrfExtentType6(inputScalars.cdrfModel6, RgEndZ, C, dz, sigmaCoeff);
        const float eEndXY = cdrfExtentType6(inputScalars.cdrfModel6, RgEndXY, C, dy, sigmaCoeff);
        extentMax = (std::max)(extentMax, (std::max)(eEndZ, eEndXY));
    }
    const int radiusStep = 8;
    int r = static_cast<int>(std::ceil(extentMax));
    r = ((r + radiusStep - 1) / radiusStep) * radiusStep; // Keeps the filter size the same across projections
    r = (std::min)(r, (std::max)(grid.Ny, grid.Nz));
    const int fSize = 2 * r + 1;
    psf.radius = r;
    if (!bank)
        return psf;

    // Distance (mm) from the collimator surface to the center of each image plane along the detector normal
    const af::array z = (af::range(af::dim4(1, 1, nPlanes), 2) + (0.5f + static_cast<float>(psf.firstPlane)) - dist) * dx;
    // Geometric collimator resolution of the two directions
    const af::array RgZ = (af::max)(Az * z + Bz, 0.f);
    const af::array RgXY = (af::max)(Axy * z + Bxy, 0.f);
    const af::array FWHMZ = af::sqrt(RgZ * RgZ + C * C) / dz; // In voxels
    const af::array FWHMXY = af::sqrt(RgXY * RgXY + C * C) / dy;
    const af::array sigmaZ = FWHMZ * sigmaCoeff;
    // The rotation of the image blurs it by roughly one voxel in the transaxial direction, which is
    // compensated for here. The axial direction is not rotated. This bank (the ArrayFire pipeline of the CPU backend) keeps the constant
    // variance sigmaCoeff^2, the GPU kernels use the one of the view (GridType6::rotVar, CDRFParams6::comp)
    const af::array sigmaXY = af::sqrt((af::max)(FWHMXY * FWHMXY - 1.f, 1e-16f)) * sigmaCoeff;

    const af::array xx = af::range(af::dim4(fSize), 0) - static_cast<float>(r);
    const af::array xx2 = af::tile(xx * xx, 1, 1, nPlanes);
    af::array gZ, gXY;
    if (inputScalars.cdrfModel6 == 1U) {
        gZ = squareBankType6(RgZ, C, dz, sigmaCoeff, false, r, nPlanes);
        gXY = squareBankType6(RgXY, C, dy, sigmaCoeff, true, r, nPlanes);
    }
    else {
        gZ = af::exp(-xx2 / af::tile(2.f * sigmaZ * sigmaZ, fSize, 1, 1));
        gXY = af::exp(-xx2 / af::tile(2.f * sigmaXY * sigmaXY, fSize, 1, 1));
    }
    // Normalize each plane separately since the Gaussians are truncated
    gZ /= af::tile(af::sum(gZ, 0), fSize, 1, 1);
    gXY /= af::tile(af::sum(gXY, 0), fSize, 1, 1);
    // Zero the planes that the rounding of firstPlane added, they are behind the collimator surface.
    // This has to be done whenever the start was rounded down, including when it was rounded to zero
    gZ *= af::tile((z > 0.f).as(f32), fSize, 1, 1);
    psf.axial = gZ;                                 // fSize x 1 x nPlanes, convolves dimension 0
    psf.transaxial = af::reorder(gXY, 1, 0, 2);     // 1 x fSize x nPlanes, convolves dimension 1
    af::eval(psf.axial);
    af::eval(psf.transaxial);
    return psf;
}

// The PSF filter bank for projection u1, either computed on-the-fly or taken from the precomputed bank
// The precomputed bank is sampled at the voxel size of the main volume and can only be shifted by an
// integer number of voxels, so it is not applicable to the multi-resolution volumes. Those always use
// the on-the-fly response, which the setup makes sure is available
inline PSFType6 getPSFType6(const Weighting& w_vec, const scalarStruct& inputScalars, const GridType6& grid, const int64_t u1, const int ii = 0) {
    // The same goes for a grid that is not the volume itself, that is an off-centre main volume, which is
    // rotated into an axis-centered grid as well
    if (onTheFlyPSFType6(inputScalars, ii)) {
        // The fused kernels compute the weights themselves, so no filter bank is built for them
        const bool fused = type6FusedAllowed(inputScalars, ii);
        PSFType6 psf = psfHelperType6(w_vec, inputScalars, grid, ii, !fused);
        psf.fused = fused;
        return psf;
    }
    /* w_vec.distInt[u1] interpretation
        * A positive value indicates that the detector is inside the FOV
        * * After shift, set first w_vec.distInt[u1] elements through z-axis to zero
        * A negative value indicates that the detector is outside the FOV
        * * After shift, set last w_vec.distInt[u1] elements through z-axis to zero
    */
    PSFType6 psf;
    af::array PSF = translateHelperType6(w_vec.gFilter, 0, 0, w_vec.distInt[u1]);
    psf.full = PSF(af::span, af::span, af::seq(0, grid.Nx - 1));
    // A user supplied filter can be of any size, so round the radius up
    psf.radius = static_cast<int>(((std::max)(psf.full.dims(0), psf.full.dims(1)) + 1) / 2);
    return psf;
}

// Convolves the image, which has to be in the (Nz, Ny, Nx) order, with the CDRF. Only the planes that
// are in front of the collimator are convolved, the rest of them do not contribute to the projection
inline af::array applyPSFType6(const af::array& im, const PSFType6& psf, const scalarStruct& inputScalars, ProjectorClass& proj) {
    if (!psf.separable)
        return af::convolve2(im, psf.full);
    af::array visible = psf.firstPlane == 0 ? im : im(af::span, af::span, af::seq(psf.firstPlane, static_cast<int>(im.dims(2)) - 1));
#ifndef CPU
    // The custom kernel is several times faster than the ArrayFire convolution, which does not
    // handle a separate filter for every slice particularly well
    visible.eval();
    af::array tmp(visible.dims(), f32);
    af::array output(visible.dims(), f32);
    if (convolveCDRFAF(tmp, visible, psf.axial, inputScalars, proj, 0) == 0
        && convolveCDRFAF(output, tmp, psf.transaxial, inputScalars, proj, 1) == 0)
        return output;
#endif
    return af::convolve2(af::convolve2(visible, psf.axial), psf.transaxial);
}

// Restores the planes that applyPSFType6 skipped. Needed by the backprojection, which has to give an
// image of the original size, but not by the forward projection, which sums over the planes anyway
inline af::array restoreDepthType6(const af::array& im, const PSFType6& psf) {
    if (!psf.separable || psf.firstPlane == 0)
        return im;
    return af::join(2, af::constant(0.f, im.dims(0), im.dims(1), psf.firstPlane), im);
}

#ifndef CPU
// Wrappers of the fused type-6 kernels, see the launchers of ProjectorClass and the kernels in
// auxKernels.cl. The output array has to be allocated by the caller, with the dimensions given there.
// The arrays are handed to the kernels as device pointers, like in convolveCDRFAF
// The response weights and the plane radii of one view and volume (nPlanes planes, the first of them is the
// plane prm.xOff of the grid), see the cdrfWeights kernel. They are made once and used by every pass of the
// view: the forward projection and both of the backprojections (measurement and sensitivity image) of it
inline int cdrfWeightsAF(PSFType6& psf, const CDRFParams6& prm, const int nPlanes, const scalarStruct& inputScalars, ProjectorClass& proj) {
    const int nTaps = 2 * prm.radius + 1;
    // The kernel writes every element
    psf.wAx = af::array(af::dim4(nPlanes, nTaps), f32);
    psf.wTr = af::array(af::dim4(nPlanes, nTaps), f32);
    psf.rAx = af::array(af::dim4(nPlanes), s32);
    psf.rTr = af::array(af::dim4(nPlanes), s32);
    proj.d_cdrfWAx = transferAF(psf.wAx);
    proj.d_cdrfWTr = transferAF(psf.wTr);
    proj.d_cdrfRAx = transferAF(psf.rAx);
    proj.d_cdrfRTr = transferAF(psf.rTr);
    const int status = proj.cdrfWeights(inputScalars, nPlanes, prm);
    psf.wAx.unlock();
    psf.wTr.unlock();
    psf.rAx.unlock();
    psf.rTr.unlock();
    return status;
}

// The axial pass. att is the attenuation of the region, or null (see AttRefType6). With
// attMode 1 it multiplies the input (forward projection), with 2 the output (backprojection). The region is N0 x N1 x N2, the output
// holds nZOut rows of it from zStart on. The input and the output may be allocated larger than that in the first two dimensions (their
// dimensions are their strides), which is what lets the arrays that vary in size from one view to the next be rounded up
inline int cdrfAxialAF(af::array& output, const af::array& input, const AttRefType6* att, const int attMode, const int N0, const int N1, const int N2,
    const int zStart, const int nZOut, const CDRFParams6& prm, const PSFType6& psf, const scalarStruct& inputScalars, ProjectorClass& proj) {
    const bool direct = prm.globalW;
    proj.d_im = transferAF(input);
    proj.d_rhs = transferAF(output);
    if (direct) {
        proj.d_cdrfWAx = transferAF(psf.wAx);
        proj.d_cdrfRAx = transferAF(psf.rAx);
    }
    CDRFLayout6 lay;
    lay.inSx = static_cast<int>(input.dims(0));
    lay.inSy = static_cast<int>(input.dims(1));
    lay.outSx = static_cast<int>(output.dims(0));
    lay.outSy = static_cast<int>(output.dims(1));
    // The kernel does not read the attenuation with attMode 0, but it needs a valid argument
    if (att != nullptr) {
        const af::array& attRef = *att->arr;
        proj.d_cdrfAtt = transferAF(attRef);
        lay.attX0 = att->x0;
        lay.attY0 = att->y0;
        lay.attZ0 = att->z0 + att->slice * att->zStride;
        lay.attSx = static_cast<int>(att->arr->dims(0));
        lay.attSy = static_cast<int>(att->arr->dims(1));
        lay.attDx = att->dx;
    }
    else
        proj.d_cdrfAtt = proj.d_im;
    const int status = proj.cdrfAxial(inputScalars, N0, N1, N2, zStart, nZOut, att != nullptr ? attMode : 0, prm, lay);
    input.unlock();
    output.unlock();
    if (direct) {
        psf.wAx.unlock();
        psf.rAx.unlock();
    }
    if (att != nullptr)
        att->arr->unlock();
    return status;
}

// The transaxial pass of the forward projection together with the sum over the depth: rows yStart ...
// yStart + nOut - 1 of the input (x, y, z) go to the view'th (nOut x input.dims(2)) slot of out. It is added
// to what is there with accumulate
inline int cdrfTransFPAF(af::array& out, const af::array& input, const int N0, const int N1, const int N2, const int yStart, const int nOut, const int view, const bool accumulate, const CDRFParams6& prm, const PSFType6& psf,
    const scalarStruct& inputScalars, ProjectorClass& proj) {
    const bool direct = prm.globalW;
    proj.d_im = transferAF(input);
    proj.d_rhs = transferAF(out);
    if (direct) {
        proj.d_cdrfWTr = transferAF(psf.wTr);
        proj.d_cdrfRTr = transferAF(psf.rTr);
    }
    const int status = proj.cdrfTransFP(inputScalars, N0, N1, N2, static_cast<int>(input.dims(0)), static_cast<int>(input.dims(1)),
        yStart, nOut, view, accumulate ? 1 : 0, prm);
    input.unlock();
    out.unlock();
    if (direct) {
        psf.wTr.unlock();
        psf.rTr.unlock();
    }
    return status;
}

// The transaxial pass of the backprojection. The projection (trans x axial) is the view'th one of the
// projections in the array, which sits at (pY0, pZ0) of the region that is the output, and is spread over
// the depth. With resample (only when the kernel makes the weights, not with prm.globalW) the array holds the measurements
// (nRowsD x nColsD x views) instead and the kernel makes the projection on the lattice of the output
// region from the view'th one of them, like resampleFromDetectorAF does
inline int cdrfTransBPAF(af::array& output, const af::array& projection, const int N0, const int N1, const int N2, const int view, const int pY0, const int pZ0, const CDRFParams6& prm, const PSFType6& psf,
    const scalarStruct& inputScalars, ProjectorClass& proj, const MapType6* resampleMap = nullptr) {
    const bool direct = prm.globalW;
    proj.d_im = transferAF(projection);
    proj.d_rhs = transferAF(output);
    if (direct) {
        proj.d_cdrfWTr = transferAF(psf.wTr);
        proj.d_cdrfRTr = transferAF(psf.rTr);
    }
    int status;
    if (resampleMap != nullptr) {
        status = proj.cdrfTransBP(inputScalars, N0, N1, N2, static_cast<int>(output.dims(0)), static_cast<int>(output.dims(1)),
            N1, N2, 0, 0, view, prm, 1, static_cast<int>(inputScalars.nRowsD),
            static_cast<int>(inputScalars.nColsD), resampleMap->trA, resampleMap->trB, resampleMap->axA, resampleMap->axB);
    }
    else {
        status = proj.cdrfTransBP(inputScalars, N0, N1, N2, static_cast<int>(output.dims(0)), static_cast<int>(output.dims(1)),
            static_cast<int>(projection.dims(0)), static_cast<int>(projection.dims(1)), pY0, pZ0, view, prm);
    }
    projection.unlock();
    output.unlock();
    if (direct) {
        psf.wTr.unlock();
        psf.rTr.unlock();
    }
    return status;
}

// Resamples the projection on the lattice of the rotated grid (trans x axial) onto the detector and adds
// it to the view'th projection of fProj (nRowsD x nColsD x views), see the resampleToDetector kernel
inline int resampleToDetectorAF(af::array& fProj, const af::array& P, const int view, const MapType6& map, const scalarStruct& inputScalars, ProjectorClass& proj) {
    proj.d_im = transferAF(P);
    proj.d_rhs = transferAF(fProj);
    const int status = proj.resampleToDetector(inputScalars, static_cast<int>(P.dims(0)), static_cast<int>(P.dims(1)),
        static_cast<int>(inputScalars.nRowsD), static_cast<int>(inputScalars.nColsD), map.trA, map.trB, map.axA, map.axB, view);
    P.unlock();
    fProj.unlock();
    return status;
}

// The adjoint of the above: the view'th measurement (nRowsD x nColsD x views) onto the lattice of the
// rotated grid, which is the output P (trans x axial)
inline int resampleFromDetectorAF(af::array& P, const af::array& meas, const int view, const MapType6& map, const scalarStruct& inputScalars, ProjectorClass& proj) {
    proj.d_im = transferAF(meas);
    proj.d_rhs = transferAF(P);
    const int status = proj.resampleFromDetector(inputScalars, static_cast<int>(P.dims(0)), static_cast<int>(P.dims(1)),
        static_cast<int>(inputScalars.nRowsD), static_cast<int>(inputScalars.nColsD), map.trA, map.trB, map.axA, map.axB, view);
    meas.unlock();
    P.unlock();
    return status;
}
#endif

// The rotation offsets (in voxels) used by the forward projection of projection u1
// The image is first rotated such that the detector is parallel with it and then shifted such that the
// center of the detector panel is at the center of the image. With swivelling detector heads the panel
// is not, in general, centered with respect to the rotation axis, in which case the shift is non-zero
// The rotation is around the rotation axis of the gantry, so the offset of a multi-resolution volume
// from the axis is part of the shift as well
inline void offsetsFPType6(const Weighting& w_vec, const GridType6& grid, const int64_t u1, float& offX, float& offY) {
    const float panelShift = grid.panelShift;
    // The rotation applied to the image by the FP is the inverse of the swivel angle
    const float cosa = std::cos(w_vec.swivelAngles[u1] * static_cast<float>(M_PI) / 180.f);
    const float sina = std::sin(w_vec.swivelAngles[u1] * static_cast<float>(M_PI) / 180.f);
    offX = -sina * panelShift + grid.axisX;
    offY = cosa * panelShift + grid.axisY;
}

// The rotation offsets (in voxels) used by the backprojection of projection u1, i.e. the inverse of
// the transformation given by offsetsFPType6. The inverse of a rotation R followed by a shift t is the
// rotation R transposed applied to -t
inline void offsetsBPType6(const Weighting& w_vec, const GridType6& grid, const int64_t u1, float& offX, float& offY) {
    float fpX = 0.f, fpY = 0.f;
    offsetsFPType6(w_vec, grid, u1, fpX, fpY);
    const float cosa = std::cos(w_vec.swivelAngles[u1] * static_cast<float>(M_PI) / 180.f);
    const float sina = std::sin(w_vec.swivelAngles[u1] * static_cast<float>(M_PI) / 180.f);
    offX = -(cosa * fpX + sina * fpY);
    offY = -(cosa * fpY - sina * fpX);
}

// The part of the rotated grid that volume ii can reach for projection u1, which the fused pipeline limits
// the work to. The footprint of the volume in the grid is worked out from the 4 corners of its xy box (the
// box of the voxel edges, -0.5 ... N - 0.5 in voxel index coordinates, which is the box of bx/by + N * d)
// with the inverse of the transformation that the rotation kernel applies, which for the output q (centered
// on the output grid) samples the source at R * q + off, i.e. q = R^T * (source - off) with the source
// centered on the source volume. Both directions of the projector (the rotation of the forward projection
// and the rotation back of the backprojection, which is its exact inverse) use the same box, so the forward
// projection and the backprojection always skip the same views and process the same planes, which keeps
// them an adjoint pair.
//   xOff, xEnd  planes (dimension 0 of the grid, the depth) xOff ... xEnd - 1 are the ones that the volume
//               reaches and that are in front of the collimator. The planes outside of them are zero in the
//               rotated volume (the bilinear footprint reaches at most 0.71 planes past the box, the margin of
//               1 plane below and 2 above is more than that) and the ones that are behind the collimator
//               have no response
//   skip        there is nothing to do: no plane, or the volume is outside the rows of the crop region
//               (widened by one row for the bilinear footprint) and so does not reach the detector
// The rows are not narrowed to the footprint, the windows of the crop region stay as they are
typedef struct DepthWindowType6_ {
    bool skip = false;
    int xOff = 0, xEnd = 0;
    // The extent of the footprint of the volume along the rows (dimension 1) of the rotated grid, in rows of the grid (the corners of its box;
    // the rows that the bilinear footprint reaches are within 0.71 rows of them), see type6WindowsType6
    double yMin = 0., yMax = 0.;
} DepthWindowType6;

inline DepthWindowType6 depthWindowType6(const Weighting& w_vec, const scalarStruct& inputScalars, const GridType6& grid, const CropType6& crop, const PSFType6& psf,
    const int64_t u1, const int ii) {
    DepthWindowType6 dw;
    float offX = 0.f, offY = 0.f;
    offsetsFPType6(w_vec, grid, u1, offX, offY);
    // The rotation of the forward projection (cosa, sina of the rotation kernel)
    const double theta = static_cast<double>(w_vec.swivelAngles[u1]) * M_PI / 180.0;
    const double cosa = std::cos(theta);
    const double sina = std::sin(theta);
    // The kernel rotates about the exact centre of the grids, N / 2 (real division), the source box is -N / 2 ... N / 2 about it
    const int sNx = static_cast<int>(inputScalars.Nx[ii]);
    const int sNy = static_cast<int>(inputScalars.Ny[ii]);
    const double sxLo = -0.5 * static_cast<double>(sNx), sxHi = 0.5 * static_cast<double>(sNx);
    const double syLo = -0.5 * static_cast<double>(sNy), syHi = 0.5 * static_cast<double>(sNy);
    double xMin = 0., xMax = 0., yMin = 0., yMax = 0.;
    for (int cc = 0; cc < 4; cc++) {
        const double dsx = ((cc & 1) ? sxHi : sxLo) - static_cast<double>(offX);
        const double dsy = ((cc & 2) ? syHi : syLo) - static_cast<double>(offY);
        // The output coordinate q of this source point, centered on the grid, and its grid index
        const double gx = (cosa * dsx + sina * dsy) + 0.5 * static_cast<double>(grid.Nx) - 0.5;
        const double gy = (-sina * dsx + cosa * dsy) + 0.5 * static_cast<double>(grid.Ny) - 0.5;
        if (cc == 0) {
            xMin = xMax = gx;
            yMin = yMax = gy;
        }
        else {
            xMin = (std::min)(xMin, gx);
            xMax = (std::max)(xMax, gx);
            yMin = (std::min)(yMin, gy);
            yMax = (std::max)(yMax, gy);
        }
    }
    dw.xOff = (std::max)(psf.firstPlane, static_cast<int>(std::floor(xMin)) - 1);
    dw.xEnd = (std::min)(grid.Nx, static_cast<int>(std::ceil(xMax)) + 2);
    const double transEnd = static_cast<double>(crop.transStart + crop.transLength);
    dw.yMin = yMin;
    dw.yMax = yMax;
    dw.skip = dw.xOff >= dw.xEnd || yMax + 1.0 < static_cast<double>(crop.transStart) || yMin - 1.0 >= transEnd;
    return dw;
}

// Builds the attenuation source of the rotation-based projector for volume ii from the host-provided
// attenuation map. The map is expected to cover the union bounding box of every multi-resolution
// volume (0..nMultiVolumes), at volume 0's own (finest) voxel size -- i.e. exactly the same "single
// dense array covering the full, pre-split FOV" convention already used for options.x0 /
// options.referenceImage with multi-resolution volumes. With nMultiVolumes == 0 the map IS the main
// volume and is used as it is, a plain reinterpretation of the host buffer. CTAttenuation must be
// true: type 6's rotate-and-integrate attenuation model needs a full image-domain map, unlike the
// precomputed LOR/sinogram-domain attenuation factors used elsewhere when CTAttenuation is false,
// which this projector cannot consume.
//
// With multi-resolution the map is resampled onto a lattice that is centered on the rotation axis and
// has the voxel size and the axial (z) lattice of volume ii, and that is large enough to hold the
// union bounding box at any angle (see attGridType6). It is that lattice that is rotated for every
// projection, which lets the attenuation be accumulated along the whole path through the material even
// though each volume is rotated and attenuated on its own
// The lattice of the source is that of volume ii, but with the axial lattice (first slice at bzV, nzV slices) given, which
// is the lattice of volume ii itself in attenuationSourceType6 and the union of the axial lattices of the volumes of an
// attenuation group in the projector that works on every volume of a view together (see attGroupsType6)
inline af::array attenuationSourceLatticeType6(const float* atten, const scalarStruct& inputScalars, const int ii, const float bzV, const int nzV) {
	float bxMap = inputScalars.bx[0], byMap = inputScalars.by[0], bzMap = inputScalars.bz[0];
	float bxHi = inputScalars.bx[0] + inputScalars.Nx[0] * inputScalars.dx[0];
	float byHi = inputScalars.by[0] + inputScalars.Ny[0] * inputScalars.dy[0];
	float bzHi = inputScalars.bz[0] + inputScalars.Nz[0] * inputScalars.dz[0];
	for (int jj = 1; jj <= inputScalars.nMultiVolumes; jj++) {
		bxMap = (std::min)(bxMap, inputScalars.bx[jj]);
		byMap = (std::min)(byMap, inputScalars.by[jj]);
		bzMap = (std::min)(bzMap, inputScalars.bz[jj]);
		bxHi = (std::max)(bxHi, inputScalars.bx[jj] + inputScalars.Nx[jj] * inputScalars.dx[jj]);
		byHi = (std::max)(byHi, inputScalars.by[jj] + inputScalars.Ny[jj] * inputScalars.dy[jj]);
		bzHi = (std::max)(bzHi, inputScalars.bz[jj] + inputScalars.Nz[jj] * inputScalars.dz[jj]);
	}
	const int NxMap = static_cast<int>(std::lround((bxHi - bxMap) / inputScalars.dx[0]));
	const int NyMap = static_cast<int>(std::lround((byHi - byMap) / inputScalars.dy[0]));
	const int NzMap = static_cast<int>(std::lround((bzHi - bzMap) / inputScalars.dz[0]));

	int gridNx = 0, gridNy = 0;
	rotGridType6(inputScalars, ii, gridNx, gridNy);
	int attNx = 0, srcNx = 0, srcNy = 0;
	attGridType6(inputScalars, ii, gridNx, gridNy, attNx, srcNx, srcNy);
	float cx = 0.f, cy = 0.f, cz = 0.f;
	axisType6(inputScalars, cx, cy, cz);
	const float bxV = cx - .5f * static_cast<float>(srcNx) * inputScalars.dx[ii], byV = cy - .5f * static_cast<float>(srcNy) * inputScalars.dy[ii];
	const bool timing = DEBUG || inputScalars.verbose >= 3;
	const auto t0 = std::chrono::steady_clock::now();
	// Summed-area table and box averages in double on the host (see resampleAttenuationType6Host), the float result uploaded once
	const std::vector<float> res = resampleAttenuationType6Host(atten, bxMap, byMap, bzMap, inputScalars.dx[0], inputScalars.dy[0], inputScalars.dz[0],
		NxMap, NyMap, NzMap, bxV, byV, bzV, inputScalars.dx[ii], inputScalars.dy[ii], inputScalars.dz[ii], srcNx, srcNy, nzV, false);
	af::array out(srcNx, srcNy, nzV, res.data());
	if (timing) {
		out.eval();
		af::sync();
		const std::chrono::duration<double> dt = std::chrono::steady_clock::now() - t0;
		mexPrintBase("Type 6 attenuation resampling of volume %d (host, double): %.3f s\n", ii, dt.count());
	}
	return out;
}

inline af::array attenuationSourceType6(const float* atten, const scalarStruct& inputScalars, const int ii) {
	if (inputScalars.nMultiVolumes == 0 && mainCenteredXYType6(inputScalars))
		return af::array(inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], atten);
	return attenuationSourceLatticeType6(atten, inputScalars, ii, inputScalars.bz[ii], static_cast<int>(inputScalars.Nz[ii]));
}

// attenuationSourceType6 of volume ii, kept in the cache of the projector while the map and the lattice are the same, instead of being built
// from the host map at every call. With the CPU backend there is no cache: the source is built into local, which is returned
// The map of the frame timestep is the one used (see attenuationFrameType6)
inline const af::array& attenuationSourceCachedType6(af::array& local, ProjectorClass& proj, const float* atten, const scalarStruct& inputScalars, const int ii,
	const int timestep = 0) {
	const AttFrameType6 map = attenuationFrameType6(inputScalars, atten, timestep);
#ifdef CPU
	// No cache with the CPU backend
	local = attenuationSourceType6(map.atten, inputScalars, ii);
	return local;
#else
	AttCacheType6& cache = attCacheType6(proj);
	const size_t slot = static_cast<size_t>(map.frame) * (static_cast<size_t>(inputScalars.nMultiVolumes) + 1ULL) + static_cast<size_t>(ii);
	if (cache.volSrc.size() <= slot) {
		cache.volSrc.resize(slot + 1);
		cache.volKey.resize(slot + 1);
	}
	const std::vector<double> key = attVolumeKeyType6(map.atten, inputScalars, ii);
	if (cache.volSrc[slot].isempty() || cache.volKey[slot] != key) {
		cache.volSrc[slot] = attenuationSourceType6(map.atten, inputScalars, ii);
		cache.volKey[slot] = key;
	}
	return cache.volSrc[slot];
#endif
}

// The attenuation factors of the voxels of a rotated attenuation map (the attenuation in 1/mm, the planes along dimension 0 from the
// detector side, voxel size dx): the accumulation along the depth and the conversion of the accumulated attenuation S_k to the factor
// as type6AttConv says (TYPE6_ATTCONV, the same as the fused kernels, see cdrfAttAt): 0 = exp(-dx * S_k), 1 = exp(-dx * (S_(k-1) + S_k) / 2),
// 2 = exp(-dx * S_(k-1)) * (1 - exp(-t)) / t, t = dx * (S_k - S_(k-1)). S_(-1) = 0
inline af::array attenuationFactorType6(const af::array& attenuation, const float dx) {
    const af::array acc = af::accum(attenuation, 0);
    const int mode = type6AttConv();
    if (mode == 0)
        return af::exp(-dx * acc);
    // The accumulated attenuation of the previous plane, zero for the first one
    const int n0 = static_cast<int>(acc.dims(0));
    af::array prev = af::constant(0.f, acc.dims(), acc.type());
    if (n0 > 1)
        prev = af::join(0, af::constant(0.f, 1, acc.dims(1), acc.dims(2), acc.dims(3), acc.type()), acc(af::seq(0, n0 - 2), af::span, af::span, af::span));
    if (mode == 1)
        return af::exp(-dx * (0.5f * (prev + acc)));
    const af::array t = dx * (acc - prev);
    // The series of g(t) = (1 - exp(-t)) / t for the small t, where 1 - exp(-t) loses the digits
    const af::array g = af::select(af::abs(t) < 0.1f, 1.f - t * (0.5f - t * (1.f / 6.f - t * (1.f / 24.f - t * (1.f / 120.f)))), (1.f - af::exp(-t)) / t);
    return af::exp(-dx * prev) * g;
}

// The attenuation of every image plane of volume ii for projection u1. The attenuation source is
// rotated and shifted like the emission volume (around the axis, since the source is centered on it),
// accumulated toward the detector and exponentiated. The result is for the same window (crop) of the
// grid as the emission volume, in the (x, y, z) order of the volume, with grid.Nx planes.
//
// The source is attNx planes wide, attNx >= grid.Nx, so that the path through all of the material is
// accumulated, see attGridType6. Plane k of the emission grid is plane k + xoff of the attenuation
// grid, xoff = (attNx - grid.Nx) / 2, since both are centered on the axis. The plane index increases
// away from the collimator, so the accumulation along dimension 0 starts at the detector side
// upload is false when the attenuation source is already in the rotation image, which it stays for all
// projections of a call
// With xEnd > xOff (the depth window of the fused pipeline) the result holds only the planes xOff ... xEnd - 1
// of the emission grid. The accumulation has to start at the detector side, so the planes from the first one of
// the attenuation grid up to the last one that is needed are rotated, and only the ones asked for are kept
// With acc and accRef (the fused pipeline, a depth window) nothing is exponentiated or cut out: the accumulated attenuation of the
// rotated window goes to *acc, which has to stay where it is as long as accRef is used, and accRef tells the axial kernel where the
// window of the emission grid is in it and the voxel size to exponentiate with, so that it reads it in place. An empty array is returned
inline af::array attenuationPlanesType6(const af::array& attenSrc, const Weighting& w_vec, const scalarStruct& inputScalars, ProjectorClass& proj,
    const GridType6& grid, const CropType6& crop, const int64_t u1, const int ii, const bool upload, const int xOffE = 0, const int xEndE = 0,
    af::array* acc = nullptr, AttRefType6* accRef = nullptr) {
    int attNx = 0, srcNx = 0, srcNy = 0;
    attGridType6(inputScalars, ii, grid.Nx, grid.Ny, attNx, srcNx, srcNy);
    const bool depthWindow = xEndE > xOffE;
    RotWindow win;
    win.srcNx = srcNx; win.srcNyFull = srcNy; win.srcNzFull = grid.Nz;
    win.outNx = attNx; win.outNyFull = grid.Ny;
    if (depthWindow) {
        win.outX = 0;
        win.outNxW = (attNx - grid.Nx) / 2 + xEndE;
    }
    win.outY = crop.transStart; win.outNy = crop.transLength;
    win.outZ = crop.axialStart; win.outNz = crop.axialLength;
    win.srcY = 0; win.srcNy = srcNy;
    // The source has the z lattice of the volume, which is placed padLo slices into the padded grid
    win.srcZ = padZLoType6(inputScalars, ii); win.srcNz = static_cast<int>(inputScalars.Nz[ii]);
    // The source is centered on the rotation axis, so it has no offset from it
    GridType6 attGrid = grid;
    attGrid.axisX = 0.f;
    attGrid.axisY = 0.f;
    float offX = 0.f, offY = 0.f;
    offsetsFPType6(w_vec, attGrid, u1, offX, offY);
    af::array attenuationImage = rotateHelperType6(attenSrc, inputScalars, proj, -w_vec.swivelAngles[u1], offX, offY, ii, 1, upload, win);
    const int xoff = (attNx - grid.Nx) / 2;
    // Accumulate the attenuation, which the fused axial kernel reads in place and exponentiates itself (see AttRefType6): the array is
    // the window that was rotated, and plane xOffE of the emission grid is plane xoff + xOffE of it
    if (accRef != nullptr && acc != nullptr) {
        *acc = af::accum(attenuationImage, 0);
        accRef->arr = acc;
        accRef->x0 = xoff + xOffE;
        accRef->y0 = 0;
        accRef->z0 = 0;
        accRef->dx = inputScalars.dx[ii];
        return af::array();
    }
    // Accumulate the attenuation, multiply each accumulated element by the voxel size and exponentiate
    attenuationImage = attenuationFactorType6(attenuationImage, inputScalars.dx[ii]);
    if (depthWindow) {
        if (xoff + xOffE != 0 || xoff + xEndE != static_cast<int>(attenuationImage.dims(0)))
            attenuationImage = attenuationImage(af::seq(xoff + xOffE, xoff + xEndE - 1), af::span, af::span);
    }
    else if (attNx != grid.Nx)
        attenuationImage = attenuationImage(af::seq(xoff, xoff + grid.Nx - 1), af::span, af::span);
    return attenuationImage;
}

// The normalization of a single projection of volume ii, the same for the forward and the back
// projection. The projection is the sum of the voxels along the ray times the length of one of them, so
// a volume with larger voxels contributes proportionally more per voxel, and the division by the area
// of a voxel of the main volume brings the scale to that of the ray-driven projectors
inline float normalizationType6(const scalarStruct& inputScalars, const int ii) {
    return (inputScalars.dx[ii] / inputScalars.dx[0]) / (inputScalars.dy[0] * inputScalars.dz[0]);
}

#ifndef CPU
// Whether the fused kernels read the response weights from global memory (made by cdrfWeights once per view and
// volume) instead of making them in shared memory: always with a radius that the latter cannot handle, and otherwise
// when the processed region (the planes of the depth window times the rows and the axial rows of the crop region) has at
// least type6WeightsThreshold voxels. It depends on the geometry of the view and volume only, so the forward
// projection, and the backprojections of the measurement and of the sensitivity image, which share the weights,
// make the same choice
inline bool type6GlobalWeights(const int radius, const DepthWindowType6& dw, const CropType6& crop) {
    if (radius > ProjectorClass::cdrfRMax)
        return true;
    const int64_t voxels = static_cast<int64_t>(dw.xEnd - dw.xOff) * static_cast<int64_t>(crop.transLength) * static_cast<int64_t>(crop.axialLength);
    return voxels >= type6WeightsThreshold;
}

// The scalars of the fused kernels for volume ii, which are the values psfHelperType6 works with. The
// planes are those of the grid, the depth window dw (and the crop region) is what the arrays hold
inline CDRFParams6 cdrfParamsType6(const scalarStruct& inputScalars, const GridType6& grid, const PSFType6& psf, const int ii, const DepthWindowType6& dw,
    const CropType6& crop) {
    CDRFParams6 prm;
    prm.xOff = dw.xOff;
    prm.globalW = type6GlobalWeights(psf.radius, dw, crop);
    prm.firstPlane = psf.firstPlane;
    prm.radius = psf.radius;
    prm.dist = grid.dist;
    prm.dxv = inputScalars.dx[ii];
    prm.dyv = inputScalars.dy[ii];
    prm.dzv = inputScalars.dz[ii];
    prm.Az = inputScalars.CDRFCoeffs6[0];
    prm.Bz = inputScalars.CDRFCoeffs6[1];
    prm.Axy = inputScalars.CDRFCoeffs6[2];
    prm.Bxy = inputScalars.CDRFCoeffs6[3];
    prm.C = inputScalars.CDRFCoeffs6[4];
    prm.sigmaCoeff = 1.f / (2.f * std::sqrt(2.f * std::log(2.f)));
    prm.comp = grid.rotVar;
    prm.norm = normalizationType6(inputScalars, ii);
    return prm;
}

// Makes the response weights of the fused kernels for one view and volume, for the planes of the depth
// window dw, when the kernels read them from global memory (type6GlobalWeights: a large region, or a radius that
// is too large for the kernels to make them themselves, ProjectorClass::cdrfRMax). Otherwise the kernels make
// them in shared memory and nothing is done here, nor for the other pipelines
inline void prepareFusedType6(PSFType6& PSF, const GridType6& grid, const scalarStruct& inputScalars, ProjectorClass& proj, const int ii, const DepthWindowType6& dw,
    const CropType6& crop) {
    if (!PSF.fused)
        return;
    const CDRFParams6 prm = cdrfParamsType6(inputScalars, grid, PSF, ii, dw, crop);
    if (!prm.globalW)
        return;
    if (cdrfWeightsAF(PSF, prm, dw.xEnd - dw.xOff, inputScalars, proj) != 0)
        mexPrint("Failed to compute the collimator response weights");
}

// Fused forward projection of one view: from the rotated (and, if att is not null, to be attenuated) region
// of the grid (x = depth, y = trans, z = axial) to the projection (trans x axial) on the detector, which is
// added to the view'th projection of fProj (nRowsD x nColsD x views). This is steps 3 to 6 of the
// ArrayFire pipeline, without the reorders. Both convolutions, the normalization, the sum over the
// depth and the resampling onto the detector are done by the kernels, and the result goes directly into
// fProj. The weights of PSF (prepareFusedType6), when there are any, have to be there
inline void forwardViewFusedType6(af::array& fProj, const int view, const af::array& kuvaRot, const AttRefType6* att, const GridType6& grid, const PSFType6& PSF, const CropType6& crop,
    const MapType6& map, const scalarStruct& inputScalars, ProjectorClass& proj, const int ii, const DepthWindowType6& dw, const int allocQ = 0) {
    // kuvaRot holds the planes of the depth window only. It may be allocated larger than that in the first two dimensions (see
    // RotWindow::allocQ), the sizes of the region are those of the window. With allocQ the array that the axial pass makes is
    // allocated with its dimensions rounded up to a multiple of it, in the same way
    const CDRFParams6 prm = cdrfParamsType6(inputScalars, grid, PSF, ii, dw, crop);
    const int N0 = dw.xEnd - dw.xOff;
    const int N1 = crop.transLength;
    const int N2 = crop.axialLength;
    if (kuvaRot.dims(0) < N0 || kuvaRot.dims(1) < N1 || kuvaRot.dims(2) < N2) {
        mexPrint("The rotated region does not match its window");
        return;
    }
    // Without multi-resolution the output is cut to the detector window here, instead of afterwards
    const bool windowed = crop.cropped && !grid.resample;
    const int zStart = windowed ? crop.axialOffset : 0;
    const int nZ = windowed ? static_cast<int>(inputScalars.nColsD) : N2;
    const int yStart = windowed ? crop.transOffset : 0;
    const int nY = windowed ? static_cast<int>(inputScalars.nRowsD) : N1;
    af::array axial(af::dim4(roundUpType6(N0, allocQ), roundUpType6(N1, allocQ), roundUpType6(nZ, allocQ)), f32);
    if (cdrfAxialAF(axial, kuvaRot, att, 1, N0, N1, N2, zStart, nZ, prm, PSF, inputScalars, proj) != 0) {
        mexPrint("Failed to apply the fused collimator response");
        return;
    }
    // The transaxial pass sums over the depth too (the planes behind the collimator are zero and are left out)
    if (!grid.resample) {
        // The projection is on the lattice of the detector, which is the size of the window
        if (nY != static_cast<int>(inputScalars.nRowsD) || nZ != static_cast<int>(inputScalars.nColsD)) {
            mexPrint("The projection does not match the detector");
            return;
        }
        if (cdrfTransFPAF(fProj, axial, N0, N1, nZ, yStart, nY, view, true, prm, PSF, inputScalars, proj) != 0)
            mexPrint("Failed to apply the fused collimator response");
        return;
    }
    // The projection (trans x axial) is resampled onto the detector: D = Wtr * P * Wax^T, see
    // toDetectorType6, and added to fProj by the resampling kernel
    af::array summed(af::dim4(nY, nZ), f32);
    if (cdrfTransFPAF(summed, axial, N0, N1, nZ, yStart, nY, 0, false, prm, PSF, inputScalars, proj) != 0
        || resampleToDetectorAF(fProj, summed, view, map, inputScalars, proj) != 0)
        mexPrint("Failed to resample the projection onto the detector");
}

// Fused backprojection of one view: from the view'th projection (trans x axial) of fp (nRowsD x nColsD x
// views) to the region of the grid, ready to be rotated back. This is steps 1 and 2 of the ArrayFire
// pipeline and the attenuation. The transaxial pass spreads the projection over the depth and normalizes
// it, the axial one convolves and attenuates. The weights of PSF (prepareFusedType6), when there are any, have to be there
inline af::array backViewFusedType6(const af::array& fp, const int view, const AttRefType6* att, const GridType6& grid, const PSFType6& PSF, const CropType6& crop,
    const MapType6& map, const scalarStruct& inputScalars, ProjectorClass& proj, const int ii, const DepthWindowType6& dw, const int allocQ = 0) {
    const CDRFParams6 prm = cdrfParamsType6(inputScalars, grid, PSF, ii, dw, crop);
    // Only the planes of the depth window, the others are not rotated back anyway
    const int N0 = dw.xEnd - dw.xOff;
    const int N1 = crop.transLength;
    const int N2 = crop.axialLength;
    const af::array* P = &fp;
    af::array resampled;
    int pView = view;
    int pY0 = 0, pZ0 = 0;
    // The measurements are resampled onto the lattice of the rotated grid inside the transaxial kernel, when it
    // can (when it makes the weights itself, see type6GlobalWeights), see the cdrfTransBP kernel
    const bool resampleInKernel = grid.resample && !prm.globalW;
    if (grid.resample && !resampleInKernel) {
        // Wtr^T * M * Wax is the transpose of fromDetectorType6
        resampled = af::array(af::dim4(N1, N2), f32);
        if (resampleFromDetectorAF(resampled, fp, view, map, inputScalars, proj) != 0)
            mexPrint("Failed to resample the measurements onto the grid");
        P = &resampled;
        pView = 0;
    }
    else if (crop.cropped) {
        // The measurements only cover the detector, the rest of the volume gets nothing
        pY0 = crop.transOffset;
        pZ0 = crop.axialOffset;
    }
    // The array that the transaxial pass makes is only read by the axial one, so it can be allocated larger (rounded up with allocQ), the
    // result goes to the rotation and is exactly N0 x N1 x N2
    af::array spread(af::dim4(roundUpType6(N0, allocQ), roundUpType6(N1, allocQ), roundUpType6(N2, allocQ)), f32);
    af::array output(af::dim4(N0, N1, N2), f32);
    if (cdrfTransBPAF(spread, *P, N0, N1, N2, pView, pY0, pZ0, prm, PSF, inputScalars, proj, resampleInKernel ? &map : nullptr) != 0
        || cdrfAxialAF(output, spread, att, 2, N0, N1, N2, 0, N2, prm, PSF, inputScalars, proj) != 0)
        mexPrint("Failed to apply the fused collimator response");
    return output;
}
#endif

// The multiplicative in-reconstruction corrections of the SPECT projector type 6: the normalization n, the multiplicative scatter /
// extra correction s (corrVector) and the global correction factor g. The projector types 1-3 multiply every ray by them (local_norm,
// local_scat, global_factor) in both the forward and the backprojection, so the model is y ~ n*s*g*(A x). The type 6 projector has no
// such kernel, so w = n .* s .* g is built here for the views uu ... uu + length - 1 (nRowsD x nColsD x length, the layout of the
// measurements of the subset) from the host arrays that the device buffers d_norm / d_scat are filled from (see
// ProjectorClass::createBuffers). The forward projection is multiplied by w, and the backprojection by w too (data is multiplied by w,
// and the sensitivity image backprojects w instead of ones), which keeps the two adjoint. Returns an empty array when no correction is
// enabled
inline af::array type6Weights(const scalarStruct& inputScalars, const Weighting& w_vec, const int64_t uu, const int64_t length, const int timestep) {
	const bool useNorm = inputScalars.normalization_correction && inputScalars.size_norm > 1ULL && w_vec.normHost != nullptr;
	const bool useScat = inputScalars.scatter == 1U && inputScalars.size_scat > 1ULL && w_vec.corrHost != nullptr;
	const bool useGlobal = inputScalars.global_factor != 1.f;
	if (!useNorm && !useScat && !useGlobal)
		return af::array();
	const size_t nDet = static_cast<size_t>(inputScalars.nRowsD) * static_cast<size_t>(inputScalars.nColsD);
	std::vector<float> wh(nDet * static_cast<size_t>(length), inputScalars.global_factor);
	if (useNorm) {
		if (inputScalars.SPECT && inputScalars.normZ == inputScalars.nHeads) {
			// One normalization image per detector head: the head of each view (the same indexing as in the kernels of the projector types 1-3)
			for (int64_t kk = 0; kk < length; kk++) {
				const size_t head = w_vec.detectorVector != nullptr ? static_cast<size_t>(w_vec.detectorVector[geomViewStartType6(w_vec, timestep) + uu + kk]) : 0ULL;
				const float* src = &w_vec.normHost[head * nDet];
				float* dst = &wh[static_cast<size_t>(kk) * nDet];
				for (size_t dd = 0; dd < nDet; dd++)
					dst[dd] *= src[dd];
			}
		}
		else {
			const float* src = &w_vec.normHost[static_cast<size_t>(normViewStartType6(w_vec, timestep) + uu) * nDet];
			for (size_t dd = 0; dd < wh.size(); dd++)
				wh[dd] *= src[dd];
		}
	}
	if (useScat) {
		// The same per-timestep frame offset as when the d_scat buffers are filled
		const size_t scatTimestepOffset = (inputScalars.size_scat > static_cast<size_t>(inputScalars.kokoNonTOF)) ?
			(static_cast<size_t>(inputScalars.kokoNonTOF) * static_cast<size_t>(timestep)) : 0ULL;
		const float* src = &w_vec.corrHost[static_cast<size_t>(uu) * nDet + scatTimestepOffset];
		for (size_t dd = 0; dd < wh.size(); dd++)
			wh[dd] *= src[dd];
	}
	return af::array(static_cast<dim_t>(inputScalars.nRowsD), static_cast<dim_t>(inputScalars.nColsD), static_cast<dim_t>(length), wh.data());
}

#ifndef CPU
// The projectors that process several views per launch (TYPE6_BATCH, see type6BatchSize), defined below with the other projectors of
// every volume. The single-volume projectors without multi-resolution go through them too
inline void forwardProjectionType6AllBatched(af::array& fProj, const Weighting& w_vec, AF_im_vectors& vec, const scalarStruct& inputScalars,
	const int64_t length, const int64_t uu, ProjectorClass& proj, const float* atten, const int timestep);
inline void backprojectionType6AllBatched(af::array& fProj, const Weighting& w_vec, AF_im_vectors& vec,
	const scalarStruct& inputScalars, const int64_t length, const int64_t uu, ProjectorClass& proj, uint32_t timestep, const uint32_t osa_iter,
	const uint32_t iter, const uint8_t compute_norm_matrix, const uint32_t iter0, const float* atten, const uint32_t sensIdx, const bool allowSensCompute);
#endif

// SPECT forward projection (projector type 6), without the multiplicative corrections (see forwardProjectionType6 below)
inline void forwardProjectionType6Raw(af::array& fProj, const Weighting& w_vec, AF_im_vectors& vec, const scalarStruct& inputScalars,
	const int64_t length, const int64_t uu, ProjectorClass& proj, const int ii = 0, const float* atten = nullptr, const int timestep = 0) {
#ifdef CPU
	if (inputScalars.nMultiVolumes > 0) {
		mexPrint("Multi-resolution volumes are not supported with the SPECT projector type 6 on the CPU backend");
		return;
	}
	if (!mainCenteredType6(inputScalars)) {
		mexPrint("An off-centre FOV is not supported with the SPECT projector type 6 on the CPU backend");
		return;
	}
#endif
#ifndef CPU
	// Several views per launch (TYPE6_BATCH): the single volume goes through the projector of every volume
	if (inputScalars.nMultiVolumes == 0 && ii == 0 && type6BatchSizeEffective(inputScalars) > 1) {
		forwardProjectionType6AllBatched(fProj, w_vec, vec, inputScalars, length, uu, proj, atten, timestep);
		return;
	}
#endif
	if (DEBUG || inputScalars.verbose >= 3) {
		af::sync();
		proj.tStartLocal = std::chrono::steady_clock::now();
	}
	if (DEBUG || inputScalars.verbose >= 3)
		mexPrint("Starting SPECT forward projection");
	int64_t u1 = uu + geomViewStartType6(w_vec, timestep);
	const af::array apuArr = ownedType6(af::moddims(vec.im_os[timestep][ii], inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii]));
	const bool useAtten = (inputScalars.attenuation_correction && inputScalars.CTAttenuation && (atten != nullptr));
	// The attenuation source does not depend on the projection, so it is built once for all of them
	af::array attenSrcLocal;
	const af::array& attenSrc = useAtten ? attenuationSourceCachedType6(attenSrcLocal, proj, atten, inputScalars, ii, static_cast<int>(timestep)) : attenSrcLocal;
	// Neither the image nor the attenuation source changes between the projections, so both are uploaded
	// only once. The first projections may be skipped entirely, so this is not simply the first one
	bool uploaded = false;

	for (int kk = 0; kk < length; kk++) {
		af::array kuvaRot, attenFused;
		AttRefType6 attFusedRef;
        const GridType6 grid = gridType6(w_vec, inputScalars, u1, ii);
		// The rotation and the translation are applied in a single interpolation
		float offX = 0.f, offY = 0.f;
		offsetsFPType6(w_vec, grid, u1, offX, offY);

        // The CDRF is needed already here since its radius determines how much of the volume can reach
        // the detector when the projections are cropped
        PSFType6 PSF = getPSFType6(w_vec, inputScalars, grid, u1, ii);
        const CropType6 crop = cropType6(inputScalars, PSF, grid, ii);
        // A multi-resolution volume that is entirely outside the axial coverage of the detector cannot
        // contribute to this projection
        if (crop.empty) {
            u1++;
            continue;
        }
        // The fused pipeline works on the planes that the volume reaches, and skips the views that it does not
        // reach at all (the same window as in the backprojection, see depthWindowType6)
        DepthWindowType6 dw;
        if (PSF.fused) {
            dw = depthWindowType6(w_vec, inputScalars, grid, crop, PSF, u1, ii);
            if (dw.skip) {
                u1++;
                continue;
            }
        }
        const MapType6 map = mapType6(inputScalars, crop, grid, ii);

        // The rotation only writes the part of the volume that can reach the detector
        const RotWindow rotWin = PSF.fused ? rotWindowOut(crop, grid, inputScalars, ii, dw.xOff, dw.xEnd) : rotWindowOut(crop, grid, inputScalars, ii);

        // 1. Rotate and translate the image
        // The same image is rotated for every projection, so it only has to be uploaded once
        kuvaRot = rotateHelperType6(apuArr, inputScalars, proj, -w_vec.swivelAngles[u1], offX, offY, ii, 0, !uploaded, rotWin);
        if (DEBUG || inputScalars.verbose > 2)
            mexPrint("Projector 6 FP steps 1 and 2 complete");

        // 3. Process and apply attenuation image
        if (useAtten) {
            // The attenuation of every plane, using the same transformation and window as above. With the
            // fused pipeline it is applied by the convolution kernel instead of here
            if (PSF.fused)
                attenuationPlanesType6(attenSrc, w_vec, inputScalars, proj, grid, crop, u1, ii, !uploaded, dw.xOff, dw.xEnd, &attenFused, &attFusedRef);
            else
                kuvaRot *= attenuationPlanesType6(attenSrc, w_vec, inputScalars, proj, grid, crop, u1, ii, !uploaded);
        }
        uploaded = true;
        if (DEBUG || inputScalars.verbose > 2)
            mexPrint("Projector 6 FP step 3 complete");
#ifndef CPU
        if (PSF.fused) {
            // Steps 4 to 6 in a few kernel launches, which add the projection to the kk'th one of fProj
            prepareFusedType6(PSF, grid, inputScalars, proj, ii, dw, crop);
            forwardViewFusedType6(fProj, kk, kuvaRot, useAtten ? &attFusedRef : nullptr, grid, PSF, crop, map, inputScalars, proj, ii, dw);
            if (DEBUG || inputScalars.verbose > 2)
                mexPrint("Projector 6 FP steps 4 to 6 complete");
            u1++;
            continue;
        }
#endif

        // 4. Convolve with detector PSF
        kuvaRot = applyPSFType6(af::reorder(kuvaRot, 2, 1, 0), PSF, inputScalars, proj);
        // Discard the margin that was only needed to make the convolution exact. With multi-resolution
        // the region is not on the lattice of the detector, so it is resampled after the sum instead
        if (crop.cropped && !grid.resample)
            kuvaRot = kuvaRot(windowAxialSeq(crop, inputScalars), windowTransSeq(crop, inputScalars), af::span);
        if (DEBUG || inputScalars.verbose > 2)
            mexPrint("Projector 6 FP step 4 complete");

        // 5. Sum over the depth, which is the last dimension after the reordering above
        kuvaRot = af::sum(kuvaRot, 2);
        kuvaRot *= normalizationType6(inputScalars, ii);
        // 6. Resample the projection onto the detector, which it is already on without multi-resolution
        if (grid.resample)
            kuvaRot = toDetectorType6(kuvaRot, crop, map, inputScalars);
        kuvaRot = af::reorder(kuvaRot, 1, 0, 2);
        if (DEBUG || inputScalars.verbose > 2)
            mexPrint("Projector 6 FP step 5 complete");

		fProj(af::span, af::span, kk) += kuvaRot;
		u1++;
	}
	if (DEBUG || inputScalars.verbose >= 3) {
		proj.tEndLocal = std::chrono::steady_clock::now();
		const std::chrono::duration<double> tDiff = proj.tEndLocal - proj.tStartLocal;
		mexPrintBase("SPECT forward projection completed in %f seconds\n", tDiff);
		mexEval();
	}
}

// SPECT forward projection (projector type 6) of volume ii, added to fProj
// Applies the multiplicative corrections (normalization, multiplicative scatter, global factor) like the projector types 1-3, i.e. the
// model is y ~ n*s*g*(A x). The projection of the volume is computed to a zeroed array first so that w multiplies only this contribution
// (fProj may already hold the other volumes) and is applied exactly once
inline void forwardProjectionType6(af::array& fProj, const Weighting& w_vec, AF_im_vectors& vec, const scalarStruct& inputScalars,
	const int64_t length, const int64_t uu, ProjectorClass& proj, const int ii = 0, const float* atten = nullptr, const int timestep = 0) {
	const af::array wgt = type6Weights(inputScalars, w_vec, uu, length, timestep);
	if (wgt.isempty()) {
		forwardProjectionType6Raw(fProj, w_vec, vec, inputScalars, length, uu, proj, ii, atten, timestep);
		return;
	}
	af::array contrib = af::constant(0.f, fProj.dims());
	forwardProjectionType6Raw(contrib, w_vec, vec, inputScalars, length, uu, proj, ii, atten, timestep);
	fProj += af::moddims(af::flat(wgt) * af::flat(contrib), fProj.dims());
}

// Backprojects a single projection angle. The PSF, the window and the attenuation depend only on the
// angle, so the caller builds them once and shares them between the measurement and the sensitivity image
// The order is the adjoint of the forward projection: spread, convolve with the PSF, attenuate and
// rotate back. The forward projection attenuates before it convolves, but the two do not commute, so it is
// the convolution that the adjoint has to apply first
// The result is added to output, the flat image of volume ii. fp holds the projections of the views
// (nRowsD x nColsD x views) and this is the view'th one. The fused pipeline rotates the region back
// straight into output, the other one adds the rotated image to it
// dw is the depth window of the fused pipeline (not used by the other one), the planes of attenuationImage
// and of the weights of PSF are those of it
inline void backProjectionType6Angle(af::array& output, const af::array& fp, const int view, const GridType6& grid, const PSFType6& PSF, const CropType6& crop, const MapType6& map,
    const AttRefType6& attenuation, const bool useAtten, const Weighting& w_vec, const scalarStruct& inputScalars, ProjectorClass& proj, const int64_t u1, const int ii,
    const DepthWindowType6& dw, const int allocQ = 0) {
    af::array kuvaRot;
    float offX = 0.f, offY = 0.f;
    offsetsBPType6(w_vec, grid, u1, offX, offY);
#ifndef CPU
    if (PSF.fused) {
        // Steps 1 to 3 in a few kernel launches, then rotate back as below. The rotation reads the region
        // and adds it to the whole image
        kuvaRot = backViewFusedType6(fp, view, useAtten ? &attenuation : nullptr, grid, PSF, crop, map, inputScalars, proj, ii, dw, allocQ);
        rotateAccumulateType6(output, kuvaRot, inputScalars, proj, w_vec.swivelAngles[u1], offX, offY, ii, rotWindowIn(crop, grid, inputScalars, ii, dw.xOff, dw.xEnd));
        return;
    }
#endif
    const af::array proj2D = fp(af::span, af::span, view);
    kuvaRot = af::reorder(proj2D, 1, 0, 2); // Transpose the projection image
    if (grid.resample) {
        // The projection of a multi-resolution volume is not on the lattice of the detector, so the
        // measurements are resampled onto the lattice of the rotated grid instead
        kuvaRot = fromDetectorType6(kuvaRot, crop, map, inputScalars);
    }
    else if (crop.cropped) {
        // The measurements only cover the detector, the rest of the volume gets nothing
        af::array plane = af::constant(0.f, crop.axialLength, crop.transLength);
        plane(windowAxialSeq(crop, inputScalars), windowTransSeq(crop, inputScalars)) = kuvaRot;
        kuvaRot = plane;
    }

    // 1. Smear the input FP across the image volume
    kuvaRot = af::tile(kuvaRot, 1, 1, grid.Nx); // Repeat through z-axis
    // The transpose of the forward projection, which weights every voxel by its own length
    kuvaRot *= normalizationType6(inputScalars, ii);

    // 2. Convolve with detector PSF
    kuvaRot = restoreDepthType6(applyPSFType6(kuvaRot, PSF, inputScalars, proj), PSF);
    // Into the image volume layout, once, instead of transposing around each step below
    kuvaRot = af::reorder(kuvaRot, 2, 1, 0);

    // 3. Attenuation correction, for the same window of the grid as the image
    // (this pipeline is not used with the attenuation of the projector that works on every volume together, which is the accumulated
    // one, see AttRefType6: the array is the attenuation factors, for the window of the volume)
    if (useAtten)
        kuvaRot *= *attenuation.arr;

    // 4. and 5. Translate and rotate the image back, i.e. invert the transformation applied by the FP
    // The rotation reads the cropped region and writes the whole image, so the region does not have
    // to be placed into a zeroed volume first
    output += af::flat(rotateHelperType6(kuvaRot, inputScalars, proj, w_vec.swivelAngles[u1], offX, offY, ii, 0, true, rotWindowIn(crop, grid, inputScalars, ii)));
}

// Computes the rotation-based projector backprojection
// Currently this computes both the backprojection itself as well as the (optional) sensitivity image at the same time
inline af::array backProjectionType6Helper(const af::array &fProj, const Weighting& w_vec, const scalarStruct &inputScalars, ProjectorClass& proj, const int64_t length, uint32_t timestep, int64_t u1, const int ii, const float* atten = nullptr,
    const af::array* fProjSens = nullptr, af::array* outputSens = nullptr) {
    const af::array fp = ownedType6(af::moddims(fProj, inputScalars.nRowsD, inputScalars.nColsD, length));
    const bool compSens = (fProjSens != nullptr && outputSens != nullptr);
    af::array fpSens;
    if (compSens)
        fpSens = ownedType6(af::moddims(*fProjSens, inputScalars.nRowsD, inputScalars.nColsD, length));
	af::array output = af::constant(0.f, inputScalars.Nx[ii] * inputScalars.Ny[ii] * inputScalars.Nz[ii], 1);
    if (compSens)
        *outputSens = af::constant(0.f, inputScalars.Nx[ii] * inputScalars.Ny[ii] * inputScalars.Nz[ii], 1);

    const bool useAtten = (inputScalars.attenuation_correction && inputScalars.CTAttenuation && (atten != nullptr));
	// Build the attenuation source once per call, not per projection angle below
	af::array attenSrcLocal;
	const af::array& attenSrc = useAtten ? attenuationSourceCachedType6(attenSrcLocal, proj, atten, inputScalars, ii, static_cast<int>(timestep)) : attenSrcLocal;
	// The attenuation source does not change between the projections, so it is uploaded only once. The
	// first projections may be skipped entirely, so this is not simply the first one
	bool attenUploaded = false;

    for (int kk = 0; kk < length; kk++) {
        const GridType6 grid = gridType6(w_vec, inputScalars, u1, ii);
        PSFType6 PSF = getPSFType6(w_vec, inputScalars, grid, u1, ii);
        const CropType6 crop = cropType6(inputScalars, PSF, grid, ii);
        if (crop.empty) {
            u1++;
            continue;
        }
        // The same planes and the same skipped views as in the forward projection, see depthWindowType6
        DepthWindowType6 dw;
        if (PSF.fused) {
            dw = depthWindowType6(w_vec, inputScalars, grid, crop, PSF, u1, ii);
            if (dw.skip) {
                u1++;
                continue;
            }
        }
        const MapType6 map = mapType6(inputScalars, crop, grid, ii);
#ifndef CPU
        // The response weights are shared by the measurement and the sensitivity image
        prepareFusedType6(PSF, grid, inputScalars, proj, ii, dw, crop);
#endif

        af::array attenuationImage;
        AttRefType6 attRef;
        if (useAtten) {
            // Rotate the attenuation source with the same transformation the forward projection applies,
            // accumulate along the projection direction and exponentiate
            if (PSF.fused)
                attenuationPlanesType6(attenSrc, w_vec, inputScalars, proj, grid, crop, u1, ii, !attenUploaded, dw.xOff, dw.xEnd, &attenuationImage, &attRef);
            else {
                attenuationImage = attenuationPlanesType6(attenSrc, w_vec, inputScalars, proj, grid, crop, u1, ii, !attenUploaded);
                attRef.arr = &attenuationImage;
            }
            attenUploaded = true;
        }
        backProjectionType6Angle(output, fp, kk, grid, PSF, crop, map, attRef, useAtten, w_vec, inputScalars, proj, u1, ii, dw);
        if (compSens)
            backProjectionType6Angle(*outputSens, fpSens, kk, grid, PSF, crop, map, attRef, useAtten, w_vec, inputScalars, proj, u1, ii, dw);
		u1++;
    }

    // The mask of every volume (see prepareMaskBPType6)
    if (w_vec.maskBP != nullptr) {
        if (!applyMaskBP(output, compSens ? outputSens : nullptr, inputScalars, w_vec, ii, static_cast<int>(timestep)))
            return af::array();
    }

    return output;
}

// SPECT backprojection (projector type 6)
// sensIdx is the vec.Summ timestep slot to use (0 for a shared sensitivity image, or the real timestep
// when inputScalars.sensPerTimestep is set -- see reconstructionAF.h). allowSensCompute additionally
// gates the sensitivity computation itself: with a shared sensitivity image and Nt > 1, only the first
// timestep processed for a given subset (tt == t0) should (re)compute/finalize it; the caller passes
// false for every later timestep so it reuses the already-finalized image untouched instead of
// overwriting it with its own (otherwise identical) computation.
inline void backprojectionType6(af::array& fProj, const Weighting& w_vec, AF_im_vectors& vec,
	const scalarStruct& inputScalars, const int64_t length, const int64_t uu, ProjectorClass& proj, uint32_t timestep, const uint32_t osa_iter = 0, const uint32_t iter = 0,
	const uint8_t compute_norm_matrix = 0, const uint32_t iter0 = 0, const int ii = 0, const float* atten = nullptr, const uint32_t sensIdx = 0, const bool allowSensCompute = true) {
#ifdef CPU
	if (inputScalars.nMultiVolumes > 0) {
		mexPrint("Multi-resolution volumes are not supported with the SPECT projector type 6 on the CPU backend");
		return;
	}
	if (!mainCenteredType6(inputScalars)) {
		mexPrint("An off-centre FOV is not supported with the SPECT projector type 6 on the CPU backend");
		return;
	}
#endif
#ifndef CPU
	// Several views per launch (TYPE6_BATCH): the single volume goes through the projector of every volume
	if (inputScalars.nMultiVolumes == 0 && ii == 0 && type6BatchSizeEffective(inputScalars) > 1) {
		backprojectionType6AllBatched(fProj, w_vec, vec, inputScalars, length, uu, proj, timestep, osa_iter, iter, compute_norm_matrix, iter0, atten, sensIdx, allowSensCompute);
		return;
	}
#endif
	if (DEBUG || inputScalars.verbose >= 3)
		mexPrint("Starting SPECT backprojection");
	const bool compSens = allowSensCompute && ((iter == iter0 && compute_norm_matrix == 2) || compute_norm_matrix == 1);
	af::array sensProj, sensOut;
	// The adjoint of the forward projection of the model y ~ n*s*g*(A x): A^T (w .* data), and the sensitivity A^T w (see type6Weights)
	const af::array wgt = type6Weights(inputScalars, w_vec, uu, length, static_cast<int>(timestep));
	const af::array fIn = wgt.isempty() ? fProj : af::array(af::flat(fProj) * af::flat(wgt));
	if (compSens) {
		if (DEBUG || inputScalars.verbose >= 3)
			mexPrint("Computing sensitivity image");
		sensProj = wgt.isempty() ? af::constant(1.f, inputScalars.nRowsD, inputScalars.nColsD, length) : wgt;
	}

    vec.rhs_os[timestep][ii] = backProjectionType6Helper(fIn, w_vec, inputScalars, proj, length, timestep, uu + geomViewStartType6(w_vec, static_cast<int64_t>(timestep)), ii, atten,
		compSens ? &sensProj : nullptr, compSens ? &sensOut : nullptr);
	vec.rhs_os[timestep][ii](vec.rhs_os[timestep][ii] < inputScalars.epps && vec.rhs_os[timestep][ii] >= 0.f) = inputScalars.epps;

	if (compSens) {
		if (compute_norm_matrix == 2) {
			vec.Summ[sensIdx][ii][osa_iter] = sensOut;
			vec.Summ[sensIdx][ii][osa_iter](vec.Summ[sensIdx][ii][osa_iter] < inputScalars.epps) = inputScalars.epps;
		} else {
			vec.Summ[sensIdx][ii][0] = sensOut;
			vec.Summ[sensIdx][ii][0](vec.Summ[sensIdx][ii][0] < inputScalars.epps) = inputScalars.epps;
		}
		if (DEBUG || inputScalars.verbose >= 3)
			mexPrint("Sensitivity image computed");
	}
	if (w_vec.maskBP != nullptr) {
		if (!resetMaskedBP(vec.rhs_os[timestep][ii], compSens ? &vec.Summ[sensIdx][ii][compute_norm_matrix == 2 ? osa_iter : 0] : nullptr, inputScalars, w_vec, ii, static_cast<int>(timestep)))
			return;
	}
	if (DEBUG || inputScalars.verbose >= 3)
		mexPrint("SPECT backprojection complete");
}
#endif

// The SPECT projectors (projector type 6) for every volume of the multi-resolution reconstruction at once. They give the
// same result as calling forwardProjectionType6 / backprojectionType6 for the volumes 0 ... nMultiVolumes in turn, and do that
// when the fused pipeline is not used (the CPU and Metal backends, the precomputed filter bank, and
// without multi-resolution, where there is just the one volume). The fused pipeline goes through the views in the outer
// loop and the volumes in the inner one, which lets the volumes that have the same attenuation source lattice share one
// rotation and accumulation of the attenuation per view (see attGroupsType6)
#if defined(METAL) || defined(CPU)
inline void forwardProjectionType6All(af::array& fProj, const Weighting& w_vec, AF_im_vectors& vec, const scalarStruct& inputScalars,
	const int64_t length, const int64_t uu, ProjectorClass& proj, const float* atten = nullptr, const int timestep = 0) {
	for (int ii = 0; ii <= static_cast<int>(inputScalars.nMultiVolumes); ii++)
		forwardProjectionType6(fProj, w_vec, vec, inputScalars, length, uu, proj, ii, atten, timestep);
}

inline void backprojectionType6All(af::array& fProj, const Weighting& w_vec, AF_im_vectors& vec,
	const scalarStruct& inputScalars, const int64_t length, const int64_t uu, ProjectorClass& proj, uint32_t timestep, const uint32_t osa_iter = 0, const uint32_t iter = 0,
	const uint8_t compute_norm_matrix = 0, const uint32_t iter0 = 0, const float* atten = nullptr, const uint32_t sensIdx = 0, const bool allowSensCompute = true) {
	for (int ii = 0; ii <= static_cast<int>(inputScalars.nMultiVolumes); ii++)
		backprojectionType6(fProj, w_vec, vec, inputScalars, length, uu, proj, timestep, osa_iter, iter, compute_norm_matrix, iter0, ii, atten, sensIdx, allowSensCompute);
}
#else
// Everything that a view needs of one volume: the grid, the response, the window and the mapping onto the detector
// (fused pipeline). active is false when the volume does not contribute to the view
typedef struct ViewVolumeType6_ {
	bool active = false;
	GridType6 grid;
	PSFType6 PSF;
	CropType6 crop;
	DepthWindowType6 dw;
	MapType6 map;
} ViewVolumeType6;

inline ViewVolumeType6 viewVolumeType6(const Weighting& w_vec, const scalarStruct& inputScalars, const int64_t u1, const int ii) {
	ViewVolumeType6 vv;
	vv.grid = gridType6(w_vec, inputScalars, u1, ii);
	vv.PSF = getPSFType6(w_vec, inputScalars, vv.grid, u1, ii);
	vv.crop = cropType6(inputScalars, vv.PSF, vv.grid, ii);
	// A multi-resolution volume that is entirely outside the axial coverage of the detector cannot contribute to the view
	if (vv.crop.empty)
		return vv;
	// The same planes and the same skipped views in both directions, see depthWindowType6
	if (vv.PSF.fused) {
		vv.dw = depthWindowType6(w_vec, inputScalars, vv.grid, vv.crop, vv.PSF, u1, ii);
		if (vv.dw.skip)
			return vv;
	}
	vv.map = mapType6(inputScalars, vv.crop, vv.grid, ii);
	vv.active = true;
	return vv;
}

// Whether the volumes are projected together: the fused pipeline for every volume, which needs more than the one volume
inline bool type6AllFused(const scalarStruct& inputScalars) {
	if (inputScalars.nMultiVolumes == 0)
		return false;
	for (int ii = 0; ii <= static_cast<int>(inputScalars.nMultiVolumes); ii++)
		if (!type6FusedAllowed(inputScalars, ii))
			return false;
	return true;
}

// The attenuation of the volumes of attenuation group grp (see attGroupsType6) for projection u1: the attenuation source of the
// group, attSrc, is rotated and accumulated once, in the frame of the group, into attAcc (the sum along the depth, which is not
// exponentiated: the axial kernel makes the factors when it reads the elements, exp(-dx * value), see AttRefType6), and every volume
// of the group that is active (vv) gets in attRefs[ii] the offset of its window (the planes of its depth window and the rows and axial
// slices of its crop region) in attAcc, so that the kernel reads it in place. The rotated window is the union of those of the
// volumes. attAcc is allocated with its dimensions rounded up (type6AllocQuantum), as is the rotated window that it is made from,
// the window itself being the logical size. Returns false when none of the volumes is active (nothing is done then). upload is false
// when the source is already in the rotation image texSlot. attAcc has to stay where it is as long as attRefs is used
inline bool attenuationGroupType6(std::vector<AttRefType6>& attRefs, af::array& attAcc, const AttGroupType6& grp, const af::array& attSrc, const std::vector<ViewVolumeType6>& vv,
	const Weighting& w_vec, const scalarStruct& inputScalars, ProjectorClass& proj, const int64_t u1, const int texSlot, const bool upload) {
	int first = -1;
	int rotX = 0, y0 = 0, y1 = 0, z0 = 0, z1 = 0;
	for (size_t mm = 0; mm < grp.members.size(); mm++) {
		const int ii = grp.members[mm];
		if (!vv[ii].active)
			continue;
		// The planes, the rows and the axial slices of the volume in the frame of the group
		const int xoff = (grp.attNx - vv[ii].grid.Nx) / 2;
		const int rowOff = (grp.rows - vv[ii].grid.Ny) / 2;
		const int zShift = grp.zOff[ii] + grp.padLo - padZLoType6(inputScalars, ii);
		const int xB = xoff + vv[ii].dw.xEnd;
		const int yA = vv[ii].crop.transStart + rowOff, yB = yA + vv[ii].crop.transLength;
		const int zA = vv[ii].crop.axialStart + zShift, zB = zA + vv[ii].crop.axialLength;
		if (first < 0) {
			first = ii;
			rotX = xB; y0 = yA; y1 = yB; z0 = zA; z1 = zB;
		}
		else {
			rotX = (std::max)(rotX, xB);
			y0 = (std::min)(y0, yA);
			y1 = (std::max)(y1, yB);
			z0 = (std::min)(z0, zA);
			z1 = (std::max)(z1, zB);
		}
	}
	if (first < 0)
		return false;
	// The accumulation has to start at the detector side, so the rotation covers the planes from the first one of the
	// attenuation grid up to the last one that any of the volumes needs
	RotWindow win;
	win.srcNx = grp.srcNx; win.srcNyFull = grp.srcNy; win.srcNzFull = grp.padLo + grp.nz + grp.padHi;
	win.outNx = grp.attNx; win.outNyFull = grp.rows;
	win.outX = 0; win.outNxW = rotX;
	win.outY = y0; win.outNy = y1 - y0;
	win.outZ = z0; win.outNz = z1 - z0;
	win.srcY = 0; win.srcNy = grp.srcNy;
	win.srcZ = grp.padLo; win.srcNz = grp.nz;
	win.allocQ = type6AllocQuantum;
	// The source is centered on the rotation axis, so it has no offset from it. The shift of the panel is the same for
	// the volumes, since it is in voxels of the same size
	GridType6 attGrid = vv[first].grid;
	attGrid.axisX = 0.f;
	attGrid.axisY = 0.f;
	float offX = 0.f, offY = 0.f;
	offsetsFPType6(w_vec, attGrid, u1, offX, offY);
	af::array attenuationImage = rotateHelperType6(attSrc, inputScalars, proj, -w_vec.swivelAngles[u1], offX, offY, first, texSlot, upload, win);
	// Accumulate the attenuation along the depth. Both arrays have the rounded up dimensions, the elements outside of the window
	// are not initialized, which does not matter: the sum is along the first dimension and each element only depends on the
	// ones before it in its own line, so the lines and the planes of the window are not touched by them
	attAcc = af::accum(attenuationImage, 0);
	for (size_t mm = 0; mm < grp.members.size(); mm++) {
		const int ii = grp.members[mm];
		if (!vv[ii].active)
			continue;
		const int xoff = (grp.attNx - vv[ii].grid.Nx) / 2;
		const int rowOff = (grp.rows - vv[ii].grid.Ny) / 2;
		const int zShift = grp.zOff[ii] + grp.padLo - padZLoType6(inputScalars, ii);
		AttRefType6& ref = attRefs[ii];
		ref.arr = &attAcc;
		ref.x0 = xoff + vv[ii].dw.xOff;
		ref.y0 = vv[ii].crop.transStart + rowOff - y0;
		ref.z0 = vv[ii].crop.axialStart + zShift - z0;
		// The attenuation factors are exp(-dx * the accumulated attenuation), made by the kernel
		ref.dx = inputScalars.dx[first];
	}
	return true;
}

// The attenuation sources of the groups, one for each, which do not depend on the view
inline std::vector<af::array> attenuationSourcesType6(const float* atten, const scalarStruct& inputScalars, const std::vector<AttGroupType6>& groups) {
	std::vector<af::array> attenSrc(groups.size());
	for (size_t gg = 0; gg < groups.size(); gg++) {
		const AttGroupType6& grp = groups[gg];
		// A single main volume that is centered on the axis uses the map as it is, which is also its own lattice
		if (inputScalars.nMultiVolumes == 0 && mainCenteredXYType6(inputScalars))
			attenSrc[gg] = attenuationSourceType6(atten, inputScalars, grp.members[0]);
		else
			attenSrc[gg] = attenuationSourceLatticeType6(atten, inputScalars, grp.members[0], grp.bz, grp.nz);
	}
	return attenSrc;
}

// The attenuation sources of the groups for a call of the projector, and which of them are in their rotation image: when the
// attenuation is used those of the cache of the projector, which are built when the map or the lattices have changed, otherwise (the
// attenuation is not used) the empty local ones
typedef struct AttSourcesType6_ {
	const std::vector<af::array>* src = nullptr;
	std::vector<char>* uploaded = nullptr;
	std::vector<af::array> localSrc;
	std::vector<char> localUploaded;
} AttSourcesType6;

// The map of the frame timestep is the one used (see attenuationFrameType6). The sources of every map are kept in the cache; the rotation images
// hold the sources of one map, so that the sources count as not uploaded when the map changes from the previous call
inline void attenuationSourcesCachedType6(AttSourcesType6& out, ProjectorClass& proj, const float* atten, const scalarStruct& inputScalars,
	const std::vector<AttGroupType6>& groups, const bool useAtten, const int timestep = 0) {
	out.localSrc.assign(groups.size(), af::array());
	out.localUploaded.assign(groups.size(), 0);
	out.src = &out.localSrc;
	out.uploaded = &out.localUploaded;
	if (!useAtten)
		return;
	const AttFrameType6 map = attenuationFrameType6(inputScalars, atten, timestep);
	AttCacheType6& cache = attCacheType6(proj);
	const size_t fr = static_cast<size_t>(map.frame);
	if (cache.src.size() <= fr) {
		cache.src.resize(fr + 1);
		cache.srcKey.resize(fr + 1);
	}
	const std::vector<double> key = attGroupKeyType6(map.atten, inputScalars, groups);
	if (cache.src[fr].size() != groups.size() || cache.srcKey[fr] != key) {
		// Another map or lattice: the sources are rebuilt, and everything that was made from the previous ones is dropped
		cache.src[fr] = attenuationSourcesType6(map.atten, inputScalars, groups);
		cache.srcKey[fr] = key;
		cache.upFrame = -1;
		cache.batches.clear();
		cache.batchBytes = 0;
	}
	if (cache.uploaded.size() != groups.size() || cache.upFrame != map.frame) {
		cache.uploaded.assign(groups.size(), 0);
		cache.upFrame = map.frame;
	}
	cache.frame = map.frame;
	out.src = &cache.src[fr];
	out.uploaded = &cache.uploaded;
}

// The attenuation of the volumes that are active (vv) in a view, attRefs[ii] (for both of the projections): every group is rotated and
// accumulated once (attenuationGroupType6), into attAcc[gg], which has to stay where it is while attRefs is used
inline void attenuationViewType6(std::vector<AttRefType6>& attRefs, std::vector<af::array>& attAcc, const std::vector<AttGroupType6>& groups, const std::vector<af::array>& attenSrc,
	const std::vector<ViewVolumeType6>& vv, std::vector<char>& attUploaded, const Weighting& w_vec, const scalarStruct& inputScalars, ProjectorClass& proj, const int64_t u1) {
	for (size_t gg = 0; gg < groups.size(); gg++) {
		if (attenuationGroupType6(attRefs, attAcc[gg], groups[gg], attenSrc[gg], vv, w_vec, inputScalars, proj, u1, rotSlotAttenuationType6 + static_cast<int>(gg), !attUploaded[gg]))
			attUploaded[gg] = 1;
	}
}

// SPECT forward projection (projector type 6) of every volume, without the multiplicative corrections (see forwardProjectionType6All below)
inline void forwardProjectionType6AllRaw(af::array& fProj, const Weighting& w_vec, AF_im_vectors& vec, const scalarStruct& inputScalars,
	const int64_t length, const int64_t uu, ProjectorClass& proj, const float* atten = nullptr, const int timestep = 0) {
	// Several views per launch (TYPE6_BATCH)
	if (type6BatchSizeEffective(inputScalars) > 1) {
		forwardProjectionType6AllBatched(fProj, w_vec, vec, inputScalars, length, uu, proj, atten, timestep);
		return;
	}
	if (!type6AllFused(inputScalars)) {
		for (int ii = 0; ii <= static_cast<int>(inputScalars.nMultiVolumes); ii++)
			forwardProjectionType6Raw(fProj, w_vec, vec, inputScalars, length, uu, proj, ii, atten, timestep);
		return;
	}
	if (DEBUG || inputScalars.verbose >= 3) {
		af::sync();
		proj.tStartLocal = std::chrono::steady_clock::now();
	}
	if (DEBUG || inputScalars.verbose >= 3)
		mexPrint("Starting SPECT forward projection of all volumes");
	const int nV = static_cast<int>(inputScalars.nMultiVolumes) + 1;
	const bool useAtten = (inputScalars.attenuation_correction && inputScalars.CTAttenuation && (atten != nullptr));
	std::vector<af::array> apuArr(nV);
	for (int ii = 0; ii < nV; ii++)
		apuArr[ii] = ownedType6(af::moddims(vec.im_os[timestep][ii], inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii]));
	// The attenuation sources do not depend on the projection, so they are built once for all of them, one for each group
	const std::vector<AttGroupType6> groups = attGroupsType6(inputScalars);
	AttSourcesType6 attSources;
	attenuationSourcesCachedType6(attSources, proj, atten, inputScalars, groups, useAtten, static_cast<int>(timestep));
	const std::vector<af::array>& attenSrc = *attSources.src;
	// Neither the images nor the attenuation sources change between the projections, so each of them is uploaded only
	// once, to a rotation image of its own. The first projections may be skipped entirely, so this is not simply the first one
	// (the attenuation sources stay in their images from a call to the next one)
	std::vector<char>& attUploaded = *attSources.uploaded;
	std::vector<char> emUploaded(nV, 0);
	int64_t u1 = uu + geomViewStartType6(w_vec, timestep);
	for (int kk = 0; kk < length; kk++) {
		std::vector<ViewVolumeType6> vv(nV);
		for (int ii = 0; ii < nV; ii++)
			vv[ii] = viewVolumeType6(w_vec, inputScalars, u1, ii);
		// The attenuation of the volumes that are active, which the convolution kernel applies (it reads it in the arrays of the groups)
		std::vector<AttRefType6> attRefs(nV);
		std::vector<af::array> attAcc(groups.size());
		if (useAtten)
			attenuationViewType6(attRefs, attAcc, groups, attenSrc, vv, attUploaded, w_vec, inputScalars, proj, u1);
		for (int ii = 0; ii < nV; ii++) {
			if (!vv[ii].active)
				continue;
			// The rotation and the translation are applied in a single interpolation, and only to the part of the volume that can
			// reach the detector
			float offX = 0.f, offY = 0.f;
			offsetsFPType6(w_vec, vv[ii].grid, u1, offX, offY);
			RotWindow rotWin = rotWindowOut(vv[ii].crop, vv[ii].grid, inputScalars, ii, vv[ii].dw.xOff, vv[ii].dw.xEnd);
			rotWin.allocQ = type6AllocQuantum;
			af::array kuvaRot = rotateHelperType6(apuArr[ii], inputScalars, proj, -w_vec.swivelAngles[u1], offX, offY, ii, rotSlotEmissionType6 + ii, !emUploaded[ii], rotWin);
			emUploaded[ii] = 1;
			// The remaining steps in a few kernel launches, which add the projection to the view'th one of fProj
			prepareFusedType6(vv[ii].PSF, vv[ii].grid, inputScalars, proj, ii, vv[ii].dw, vv[ii].crop);
			forwardViewFusedType6(fProj, kk, kuvaRot, useAtten ? &attRefs[ii] : nullptr, vv[ii].grid, vv[ii].PSF, vv[ii].crop, vv[ii].map, inputScalars, proj, ii, vv[ii].dw,
				type6AllocQuantum);
		}
		if (DEBUG || inputScalars.verbose > 2)
			mexPrint("Projector 6 FP of all volumes, view complete");
		u1++;
	}
	if (DEBUG || inputScalars.verbose >= 3) {
		proj.tEndLocal = std::chrono::steady_clock::now();
		const std::chrono::duration<double> tDiff = proj.tEndLocal - proj.tStartLocal;
		mexPrintBase("SPECT forward projection of all volumes completed in %f seconds\n", tDiff);
		mexEval();
	}
}

// SPECT forward projection (projector type 6) of every volume, added to fProj
// Applies the multiplicative corrections like the projector types 1-3 (y ~ n*s*g*(A x)): the sum over the volumes is computed to a zeroed
// array and multiplied by w once (see type6Weights)
inline void forwardProjectionType6All(af::array& fProj, const Weighting& w_vec, AF_im_vectors& vec, const scalarStruct& inputScalars,
	const int64_t length, const int64_t uu, ProjectorClass& proj, const float* atten = nullptr, const int timestep = 0) {
	const af::array wgt = type6Weights(inputScalars, w_vec, uu, length, timestep);
	if (wgt.isempty()) {
		forwardProjectionType6AllRaw(fProj, w_vec, vec, inputScalars, length, uu, proj, atten, timestep);
		return;
	}
	af::array contrib = af::constant(0.f, fProj.dims());
	forwardProjectionType6AllRaw(contrib, w_vec, vec, inputScalars, length, uu, proj, atten, timestep);
	fProj += af::moddims(af::flat(wgt) * af::flat(contrib), fProj.dims());
}

// SPECT backprojection (projector type 6) of every volume. The arguments are those of backprojectionType6 (see there for sensIdx
// and allowSensCompute), without the volume. Every volume gets its own vec.rhs_os[timestep][ii], and, when the sensitivity image
// is computed, vec.Summ[sensIdx][ii][...]
inline void backprojectionType6All(af::array& fProj, const Weighting& w_vec, AF_im_vectors& vec,
	const scalarStruct& inputScalars, const int64_t length, const int64_t uu, ProjectorClass& proj, uint32_t timestep, const uint32_t osa_iter = 0, const uint32_t iter = 0,
	const uint8_t compute_norm_matrix = 0, const uint32_t iter0 = 0, const float* atten = nullptr, const uint32_t sensIdx = 0, const bool allowSensCompute = true) {
	// Several views per launch (TYPE6_BATCH)
	if (type6BatchSizeEffective(inputScalars) > 1) {
		backprojectionType6AllBatched(fProj, w_vec, vec, inputScalars, length, uu, proj, timestep, osa_iter, iter, compute_norm_matrix, iter0, atten, sensIdx, allowSensCompute);
		return;
	}
	if (!type6AllFused(inputScalars)) {
		for (int ii = 0; ii <= static_cast<int>(inputScalars.nMultiVolumes); ii++)
			backprojectionType6(fProj, w_vec, vec, inputScalars, length, uu, proj, timestep, osa_iter, iter, compute_norm_matrix, iter0, ii, atten, sensIdx, allowSensCompute);
		return;
	}
	if (DEBUG || inputScalars.verbose >= 3)
		mexPrint("Starting SPECT backprojection of all volumes");
	const int nV = static_cast<int>(inputScalars.nMultiVolumes) + 1;
	const bool compSens = allowSensCompute && ((iter == iter0 && compute_norm_matrix == 2) || compute_norm_matrix == 1);
	af::array sensProj;
	// The adjoint of the forward projection of the model y ~ n*s*g*(A x): A^T (w .* data), and the sensitivity A^T w (see type6Weights)
	const af::array wgt = type6Weights(inputScalars, w_vec, uu, length, static_cast<int>(timestep));
	const af::array fIn = wgt.isempty() ? fProj : af::array(af::flat(fProj) * af::flat(wgt));
	if (compSens) {
		if (DEBUG || inputScalars.verbose >= 3)
			mexPrint("Computing sensitivity image");
		sensProj = wgt.isempty() ? af::constant(1.f, inputScalars.nRowsD, inputScalars.nColsD, length) : wgt;
	}
	const af::array fp = ownedType6(af::moddims(fIn, inputScalars.nRowsD, inputScalars.nColsD, length));
	af::array fpSens;
	if (compSens)
		fpSens = ownedType6(af::moddims(sensProj, inputScalars.nRowsD, inputScalars.nColsD, length));
	std::vector<af::array> output(nV), outputSens(nV);
	for (int ii = 0; ii < nV; ii++) {
		output[ii] = af::constant(0.f, inputScalars.Nx[ii] * inputScalars.Ny[ii] * inputScalars.Nz[ii], 1);
		if (compSens)
			outputSens[ii] = af::constant(0.f, inputScalars.Nx[ii] * inputScalars.Ny[ii] * inputScalars.Nz[ii], 1);
	}
	const bool useAtten = (inputScalars.attenuation_correction && inputScalars.CTAttenuation && (atten != nullptr));
	// The attenuation sources are built once, one for each group, and uploaded only once (for the whole reconstruction)
	const std::vector<AttGroupType6> groups = attGroupsType6(inputScalars);
	AttSourcesType6 attSources;
	attenuationSourcesCachedType6(attSources, proj, atten, inputScalars, groups, useAtten, static_cast<int>(timestep));
	const std::vector<af::array>& attenSrc = *attSources.src;
	std::vector<char>& attUploaded = *attSources.uploaded;
	const AttRefType6 noAtten;
	int64_t u1 = uu + geomViewStartType6(w_vec, timestep);
	for (int kk = 0; kk < length; kk++) {
		std::vector<ViewVolumeType6> vv(nV);
		for (int ii = 0; ii < nV; ii++) {
			vv[ii] = viewVolumeType6(w_vec, inputScalars, u1, ii);
			// The response weights are shared by the measurement and the sensitivity image
			if (vv[ii].active)
				prepareFusedType6(vv[ii].PSF, vv[ii].grid, inputScalars, proj, ii, vv[ii].dw, vv[ii].crop);
		}
		std::vector<AttRefType6> attRefs(nV);
		std::vector<af::array> attAcc(groups.size());
		// The same transformation as in the forward projection, accumulated along the projection direction (the kernel exponentiates)
		if (useAtten)
			attenuationViewType6(attRefs, attAcc, groups, attenSrc, vv, attUploaded, w_vec, inputScalars, proj, u1);
		for (int ii = 0; ii < nV; ii++) {
			if (!vv[ii].active)
				continue;
			const AttRefType6& attenuation = useAtten ? attRefs[ii] : noAtten;
			backProjectionType6Angle(output[ii], fp, kk, vv[ii].grid, vv[ii].PSF, vv[ii].crop, vv[ii].map, attenuation, useAtten, w_vec, inputScalars, proj, u1, ii, vv[ii].dw,
				type6AllocQuantum);
			if (compSens)
				backProjectionType6Angle(outputSens[ii], fpSens, kk, vv[ii].grid, vv[ii].PSF, vv[ii].crop, vv[ii].map, attenuation, useAtten, w_vec, inputScalars, proj, u1, ii, vv[ii].dw,
					type6AllocQuantum);
		}
		u1++;
	}
	for (int ii = 0; ii < nV; ii++) {
		// The mask of every volume (see prepareMaskBPType6)
		if (w_vec.maskBP != nullptr) {
			if (!applyMaskBP(output[ii], compSens ? &outputSens[ii] : nullptr, inputScalars, w_vec, ii, static_cast<int>(timestep)))
				return;
		}
		vec.rhs_os[timestep][ii] = output[ii];
		vec.rhs_os[timestep][ii](vec.rhs_os[timestep][ii] < inputScalars.epps && vec.rhs_os[timestep][ii] >= 0.f) = inputScalars.epps;
		if (compSens) {
			if (compute_norm_matrix == 2) {
				vec.Summ[sensIdx][ii][osa_iter] = outputSens[ii];
				vec.Summ[sensIdx][ii][osa_iter](vec.Summ[sensIdx][ii][osa_iter] < inputScalars.epps) = inputScalars.epps;
			}
			else {
				vec.Summ[sensIdx][ii][0] = outputSens[ii];
				vec.Summ[sensIdx][ii][0](vec.Summ[sensIdx][ii][0] < inputScalars.epps) = inputScalars.epps;
			}
			if (DEBUG || inputScalars.verbose >= 3)
				mexPrint("Sensitivity image computed");
		}
		if (w_vec.maskBP != nullptr) {
			if (!resetMaskedBP(vec.rhs_os[timestep][ii], compSens ? &vec.Summ[sensIdx][ii][compute_norm_matrix == 2 ? osa_iter : 0] : nullptr, inputScalars, w_vec, ii, static_cast<int>(timestep)))
				return;
		}
	}
	if (DEBUG || inputScalars.verbose >= 3)
		mexPrint("SPECT backprojection of all volumes complete");
}

// ---------------------------------------------------------------------------------------------------------------------------------
// Several views per launch (TYPE6_BATCH, see type6BatchSize). The fused pipeline above launches the kernels once for every view
// of every volume, which is ~10 launches per view and volume. Here the views of a batch (consecutive views, B at most) are processed
// together, one volume at a time: the rotation, the axial pass, the transaxial pass (and the resampling, or the rotation back) are one
// launch each for the views of the batch that the volume reaches (the skipped views are left out). The arrays of the per-view stages are
// stacked, one slice per view, with the extents of the largest view rounded up (type6AllocQuantum), and the geometry of each view (the
// windows, the offsets, the radius and first plane, the attenuation, the slot of the projection) is in a block of floats of a parameter
// array that is uploaded once for the batch of a volume (type6BatchStride floats per view, see precomp.h). Each output element is
// computed with the arithmetic of the unbatched kernels: the forward projection of a view only writes its own slice and projection, and
// the rotation back adds the views of the batch in their order, in the thread of the voxel, to the value that the voxel has (rotateBackB),
// so there are no atomics and the result is independent of the batch size up to what the compiler does to the sums
//  - the views of a volume that use the global weights (type6GlobalWeights: a large region or a radius that the kernels cannot handle,
//    which is also where a launch has plenty of work) are not batched, they are processed one by one as without batching, in the
//    place that they have in the order of the views. The batched kernels make the weights in shared memory
//  - the rotation back always reads its source from a buffer, also with images (the sources of the views are the slices of one array,
//    not a texture), so with images (useBuffers = false) it is the bilinear interpolation of the buffer and not that of the texture
//    hardware; the forward rotation uses the texture of the volume (and of the attenuation source) as the unbatched one
//  - the attenuation of an attenuation group is rotated for the views of the batch with one launch, and accumulated with one af::accum of
//    the stack; the volumes read their slice of it (AttRefType6::slice)

// The part of the rotation image of a slot that the unbatched rotation uploads (see rotateCustomAF), without rotating. Nothing
// is done with buffers
inline int rotUploadType6(const af::array& im, RotWindow win, const scalarStruct& inputScalars, ProjectorClass& proj, const int texSlot) {
	if (inputScalars.useBuffers)
		return 0;
	int status = 0;
	if (win.outNxW <= 0)
		win.outNxW = win.outNx;
	if (win.srcNxB <= 0)
		win.srcNxB = win.srcNx;
#if defined(CUDA) || defined(HIP)
	AFDEVBUFF_t input = transferAF(im);
	status = proj.transferRotTex(inputScalars, input, win, texSlot);
	im.unlock();
	if (status != 0)
		return -1;
#elif defined(OPENCL)
	if (proj.ensureRotImage(win, texSlot) != 0)
		return -1;
	const cl::detail::size_t_array rotOrigin = { { static_cast<::size_t>(win.srcX), static_cast<::size_t>(win.srcY), static_cast<::size_t>(win.srcZ) } };
	const cl::detail::size_t_array rotRegion = { { static_cast<::size_t>(win.srcNxB), static_cast<::size_t>(win.srcNy), static_cast<::size_t>(win.srcNz) } };
	status = proj.CLCommandQueue[0].enqueueCopyBufferToImage(cl::Buffer(*im.device<cl_mem>(), true), proj.d_rotI[texSlot], 0, rotOrigin, rotRegion);
	im.unlock();
	if (status != 0) {
		getErrorString(status);
		mexPrint("Failed to copy rotation image\n");
		return -1;
	}
#endif
	return 0;
}

// The rotation part of the parameter block of a view (the values that rotateCustomAF passes to the rotate kernel). angleDeg is the
// angle that rotateHelperType6 is given (-swivel for the forward projection, +swivel for the rotation back)
inline void type6PutRotation(float* p, const RotWindow& w, const float angleDeg, const float offX, const float offY) {
	const float angleRad = angleDeg * static_cast<float>(M_PI) / 180.f;
	const int outNxW = w.outNxW > 0 ? w.outNxW : w.outNx;
	const int srcNxB = w.srcNxB > 0 ? w.srcNxB : w.srcNx;
	p[type6P::SRCNX] = static_cast<float>(w.srcNx);
	p[type6P::SRCNY] = static_cast<float>(w.srcNyFull);
	p[type6P::SRCNZ] = static_cast<float>(w.srcNzFull);
	p[type6P::OUTNX] = static_cast<float>(w.outNx);
	p[type6P::OUTNYF] = static_cast<float>(w.outNyFull);
	p[type6P::COSA] = std::cos(-angleRad);
	p[type6P::SINA] = std::sin(-angleRad);
	p[type6P::OFFX] = offX;
	p[type6P::OFFY] = offY;
	p[type6P::OUTY] = static_cast<float>(w.outY);
	p[type6P::OUTZ] = static_cast<float>(w.outZ);
	p[type6P::OUTNY] = static_cast<float>(w.outNy);
	p[type6P::OUTNZ] = static_cast<float>(w.outNz);
	p[type6P::SRCY] = static_cast<float>(w.srcY);
	p[type6P::SRCZ] = static_cast<float>(w.srcZ);
	p[type6P::SRCNYB] = static_cast<float>(w.srcNy);
	p[type6P::SRCNZB] = static_cast<float>(w.srcNz);
	p[type6P::OUTX] = static_cast<float>(w.outX);
	p[type6P::OUTNXW] = static_cast<float>(outNxW);
	p[type6P::SRCX] = static_cast<float>(w.srcX);
	p[type6P::SRCNXB] = static_cast<float>(srcNxB);
}

// One view of a volume in a batch: the geometry that the kernels need, worked out like the unbatched functions do
typedef struct BatchItemType6_ {
	// The slot of the view in the projections and the measurements, and its projection (angle) index
	int kk = 0;
	int64_t u1 = 0;
	ViewVolumeType6* vv = nullptr;
	// The attenuation of the view (null without), see AttRefType6
	const AttRefType6* att = nullptr;
	CDRFParams6 prm;
	// The region (the planes of the depth window, the rows and the axial rows of the crop region)
	int N0 = 0, N1 = 0, N2 = 0;
	// The forward projection keeps the rows yStart ... yStart + nY - 1 and the axial rows zStart ... zStart + nZ - 1 of the region
	// (the whole region, unless the detector window is cut out). The rotation back has nY = N1 and nZ = N2
	int yStart = 0, nY = 0, zStart = 0, nZ = 0;
	float offX = 0.f, offY = 0.f;
	// The window of the rotation (rotWindowOut for the forward projection, rotWindowIn for the rotation back)
	RotWindow win;
	// The windows of the batched kernels (see type6WindowsType6): the rows of the region along y that the stacked arrays hold (y0, ny), the
	// axial rows that the input of the axial pass holds (iz0, inz), and the window of the rotation of the batched launches (bwin, win is the
	// one of the unbatched pipeline). For the views that are not windowed they are the whole region and win
	int y0 = 0, ny = 0, iz0 = 0, inz = 0;
	RotWindow bwin;
	// Nothing of the view reaches the detector (or the image): it is left out, which is the same as adding zeros
	bool empty = false;
} BatchItemType6;

// The work windows of the multi-resolution volumes: the batched kernels only process the elements that take part in the result, which are
// exactly the same as those of the processing of the whole region. The output rows of the passes are the rows that the detector reads
// (type6TouchedType6), and the axial rows that the backprojection reads; the rotation and the stacked arrays are limited to the footprint
// of the volume, and the rotation and the accumulation of the attenuation to what the kernels of the volumes of the group read
// (type6AttWindowType6)

// The rows lo ... hi - 1 of a lattice of n samples that resampleToDetectorB (and its adjoint, cdrfFromDetectorAt) touches for the nDet detector
// bins at the positions A * d + B: for every bin the two samples floor(t) and floor(t) + 1, also when the weight of the second one is zero.
// The positions are made in double from the values that the kernels get, with a margin that is far larger than the rounding of the
// kernels (float, possibly fused), so this is never fewer rows than the kernel reads with a nonzero weight
inline void type6TouchedType6(const int nDet, const float A, const float B, const int n, int& lo, int& hi) {
	const double eps = 1e-3;
	const double t0 = static_cast<double>(B);
	const double t1 = static_cast<double>(A) * static_cast<double>(nDet - 1) + static_cast<double>(B);
	const double tMin = (std::min)(t0, t1);
	const double tMax = (std::max)(t0, t1);
	lo = (std::max)(static_cast<int>(std::floor(tMin - eps)), 0);
	hi = (std::min)(static_cast<int>(std::floor(tMax + eps)) + 2, n);
}

// The windows of the batched kernels of one view of a multi-resolution volume (the resampled projection, see grid.resample), which only
// differ from the whole region for the views that are windowed (see below). The forward projection and the backprojection use the same sets, which
// keeps them an adjoint pair:
//  - the projection that the detector reads (ty0 ... ty1 - 1 of the rows, tz0 ... tz1 - 1 of the axial rows of the lattice of the rotated
//    grid): the forward projection only makes those rows (the output of the axial pass and of the transaxial one), the lattice
//    projection of the backprojection is zero outside of them
//  - the axial rows that the volume has (its rows, padLo ... padLo + Nz - 1 of the grid, the others are padding): the backprojection
//    only makes the axial rows (output of the axial pass) that the rotation back reads, which are those, and that are within the radius of
//    the response from the rows of the projection (the others are zero)
// The views that are left empty by this are marked, they add nothing. The views that use the global weights are not batched and keep the
// whole region (their path does not read the windows)
// The windows themselves, in region coordinates of the volume (see type6WindowsType6): windowed is false when the view keeps the whole region
// (the unbatched views and the views of the volumes that are not resampled), empty when nothing of the view reaches the detector or the image
typedef struct WinRegionType6_ {
	bool windowed = false, empty = false, yWindow = false;
	int ty0 = 0, ty1 = 0, tz0 = 0, tz1 = 0;
	int wy0 = 0, wy1 = 0;
	int oz0 = 0, oz1 = 0;
} WinRegionType6;

inline WinRegionType6 type6WinRegionType6(const ViewVolumeType6& vv, const scalarStruct& inputScalars, const int ii, const bool globalW, const int N1, const int N2) {
	WinRegionType6 w;
	w.wy0 = 0; w.wy1 = N1; w.oz0 = 0; w.oz1 = N2;
	if (!vv.grid.resample || globalW)
		return w;
	w.windowed = true;
	const CropType6& crop = vv.crop;
	type6TouchedType6(static_cast<int>(inputScalars.nRowsD), vv.map.trA, vv.map.trB, N1, w.ty0, w.ty1);
	type6TouchedType6(static_cast<int>(inputScalars.nColsD), vv.map.axA, vv.map.axB, N2, w.tz0, w.tz1);
	const int r = vv.PSF.radius;
	const int vz0 = padZLoType6(inputScalars, ii) - crop.axialStart;
	const int vz1 = vz0 + static_cast<int>(inputScalars.Nz[ii]);
	w.oz0 = (std::max)((std::max)(w.tz0 - r, 0), vz0);
	w.oz1 = (std::min)((std::min)(w.tz1 + r, N2), vz1);
	// The rows of the stacked arrays: the footprint of the volume (the rotation is zero outside of it, with the margin of the depth window) within
	// the radius of the response from the rows of the projection. The elements outside of them are zero
	if (vv.PSF.fused) {
		w.wy0 = (std::max)((std::max)(static_cast<int>(std::floor(vv.dw.yMin)) - 1 - crop.transStart, w.ty0 - r), 0);
		w.wy1 = (std::min)((std::min)(static_cast<int>(std::ceil(vv.dw.yMax)) + 2 - crop.transStart, w.ty1 + r), N1);
		w.yWindow = true;
	}
	if (w.ty1 <= w.ty0 || w.tz1 <= w.tz0 || w.oz1 <= w.oz0 || w.wy1 <= w.wy0)
		w.empty = true;
	return w;
}

// The rows and the axial rows of the region of the volume that the batched kernels read of its attenuation (the forward projection reads the
// axial rows oz0 ... oz1 - 1 in the input of the axial pass, the backprojection the same ones in its output, both the rows wy0 ... wy1 - 1), for the
// attenuation windows: the whole region unless the view is windowed, and always the whole region of the views that are not
// batched. Returns false when the view of the volume is empty, which reads nothing (its items are left out)
inline bool type6AttWindowType6(const ViewVolumeType6& vv, const scalarStruct& inputScalars, const int ii, int& y0, int& y1, int& z0, int& z1) {
	y0 = 0; y1 = vv.crop.transLength; z0 = 0; z1 = vv.crop.axialLength;
	const CDRFParams6 prm = cdrfParamsType6(inputScalars, vv.grid, vv.PSF, ii, vv.dw, vv.crop);
	const WinRegionType6 w = type6WinRegionType6(vv, inputScalars, ii, prm.globalW, y1, z1);
	if (!w.windowed)
		return true;
	if (w.empty)
		return false;
	y0 = w.wy0; y1 = w.wy1; z0 = w.oz0; z1 = w.oz1;
	return true;
}

inline void type6WindowsType6(BatchItemType6& it, const scalarStruct& inputScalars, const int ii, const bool forward) {
	const ViewVolumeType6& vv = *it.vv;
	const WinRegionType6 w = type6WinRegionType6(vv, inputScalars, ii, it.prm.globalW, it.N1, it.N2);
	if (!w.windowed)
		return;
	const CropType6& crop = vv.crop;
	const int ty0 = w.ty0, ty1 = w.ty1, tz0 = w.tz0, tz1 = w.tz1, oz0 = w.oz0, oz1 = w.oz1;
	if (w.empty) {
		it.empty = true;
		return;
	}
	if (w.yWindow) {
		it.y0 = w.wy0;
		it.ny = w.wy1 - w.wy0;
	}
	if (forward) {
		it.yStart = ty0;
		it.nY = ty1 - ty0;
		it.zStart = tz0;
		it.nZ = tz1 - tz0;
		// The rotation makes the window and the axial pass reads it
		it.iz0 = oz0;
		it.inz = oz1 - oz0;
		it.bwin.outY = crop.transStart + it.y0;
		it.bwin.outNy = it.ny;
		it.bwin.outZ = crop.axialStart + oz0;
		it.bwin.outNz = oz1 - oz0;
	}
	else {
		it.bwin.srcY = crop.transStart + it.y0;
		it.bwin.srcNy = it.ny;
		it.zStart = oz0;
		it.nZ = oz1 - oz0;
		it.iz0 = tz0;
		it.inz = tz1 - tz0;
		it.bwin.srcZ = crop.axialStart + oz0;
		it.bwin.srcNz = oz1 - oz0;
	}
}

inline void type6ItemFP(BatchItemType6& it, ViewVolumeType6& vv, const AttRefType6* att, const int kk, const int64_t u1, const Weighting& w_vec, const scalarStruct& inputScalars,
	const int ii) {
	it.kk = kk;
	it.u1 = u1;
	it.vv = &vv;
	it.att = att;
	it.prm = cdrfParamsType6(inputScalars, vv.grid, vv.PSF, ii, vv.dw, vv.crop);
	it.N0 = vv.dw.xEnd - vv.dw.xOff;
	it.N1 = vv.crop.transLength;
	it.N2 = vv.crop.axialLength;
	// Without multi-resolution the output is cut to the detector window, see forwardViewFusedType6
	const bool windowed = vv.crop.cropped && !vv.grid.resample;
	it.zStart = windowed ? vv.crop.axialOffset : 0;
	it.nZ = windowed ? static_cast<int>(inputScalars.nColsD) : it.N2;
	it.yStart = windowed ? vv.crop.transOffset : 0;
	it.nY = windowed ? static_cast<int>(inputScalars.nRowsD) : it.N1;
	offsetsFPType6(w_vec, vv.grid, u1, it.offX, it.offY);
	it.win = rotWindowOut(vv.crop, vv.grid, inputScalars, ii, vv.dw.xOff, vv.dw.xEnd);
	it.y0 = 0; it.ny = it.N1; it.iz0 = 0; it.inz = it.N2;
	it.bwin = it.win;
	type6WindowsType6(it, inputScalars, ii, true);
}

inline void type6ItemBP(BatchItemType6& it, ViewVolumeType6& vv, const AttRefType6* att, const int kk, const int64_t u1, const Weighting& w_vec, const scalarStruct& inputScalars,
	const int ii) {
	it.kk = kk;
	it.u1 = u1;
	it.vv = &vv;
	it.att = att;
	it.prm = cdrfParamsType6(inputScalars, vv.grid, vv.PSF, ii, vv.dw, vv.crop);
	it.N0 = vv.dw.xEnd - vv.dw.xOff;
	it.N1 = vv.crop.transLength;
	it.N2 = vv.crop.axialLength;
	it.yStart = 0;
	it.nY = it.N1;
	it.zStart = 0;
	it.nZ = it.N2;
	offsetsBPType6(w_vec, vv.grid, u1, it.offX, it.offY);
	it.win = rotWindowIn(vv.crop, vv.grid, inputScalars, ii, vv.dw.xOff, vv.dw.xEnd);
	it.y0 = 0; it.ny = it.N1; it.iz0 = 0; it.inz = it.N2;
	it.bwin = it.win;
	type6WindowsType6(it, inputScalars, ii, false);
}

// The parameters that the passes share: the region, the response and the attenuation of the view, and its map onto the detector
inline void type6PutView(float* p, const BatchItemType6& it) {
	p[type6P::N0] = static_cast<float>(it.N0);
	p[type6P::N1] = static_cast<float>(it.N1);
	p[type6P::N2] = static_cast<float>(it.N2);
	p[type6P::ZSTART] = static_cast<float>(it.zStart);
	p[type6P::NZOUT] = static_cast<float>(it.nZ);
	p[type6P::XOFF] = static_cast<float>(it.prm.xOff);
	p[type6P::FIRSTPLANE] = static_cast<float>(it.prm.firstPlane);
	p[type6P::RADIUS] = static_cast<float>(it.prm.radius);
	p[type6P::DIST] = it.prm.dist;
	if (it.att != nullptr) {
		p[type6P::ATTX0] = static_cast<float>(it.att->x0);
		p[type6P::ATTY0] = static_cast<float>(it.att->y0);
		// The attenuation array is a stack of the attenuations of the views: the slice of the view starts zStride planes along z apart
		p[type6P::ATTZ0] = static_cast<float>(it.att->z0 + it.att->slice * it.att->zStride);
		p[type6P::ATTDX] = it.att->dx;
	}
	p[type6P::YSTART] = static_cast<float>(it.yStart);
	p[type6P::NOUT] = static_cast<float>(it.nY);
	p[type6P::VIEW] = static_cast<float>(it.kk);
	p[type6P::NTR] = static_cast<float>(it.nY);
	p[type6P::NAX] = static_cast<float>(it.nZ);
	p[type6P::TRA] = it.vv->map.trA;
	p[type6P::TRB] = it.vv->map.trB;
	p[type6P::AXA] = it.vv->map.axA;
	p[type6P::AXB] = it.vv->map.axB;
	p[type6P::Y0] = static_cast<float>(it.y0);
	p[type6P::NY] = static_cast<float>(it.ny);
	p[type6P::IZ0] = static_cast<float>(it.iz0);
	p[type6P::INZ] = static_cast<float>(it.inz);
	p[type6P::COMP] = it.prm.comp;
}

// The extents of the views of a run: the largest of each of the sizes of the regions (the largest number of the rows and of the
// axial rows of the projections too)
typedef struct RunExtType6_ {
	int n0 = 0, n1 = 0, n2 = 0, ny = 0, nz = 0;
	// The largest window along y of the stacked arrays and the largest input window of the axial pass (see BatchItemType6)
	int wy = 0, iz = 0;
} RunExtType6;

inline void type6RunInclude(RunExtType6& e, const BatchItemType6& it) {
	e.n0 = (std::max)(e.n0, it.N0);
	e.n1 = (std::max)(e.n1, it.N1);
	e.n2 = (std::max)(e.n2, it.N2);
	e.ny = (std::max)(e.ny, it.nY);
	e.nz = (std::max)(e.nz, it.nZ);
	e.wy = (std::max)(e.wy, it.ny);
	e.iz = (std::max)(e.iz, it.inz);
}

// The parameter blocks of all the launches of a batch of views (the attenuation of its groups, the runs of its volumes) are planned into one
// host vector first and uploaded with one copy (see beginBatchParams in ProjectorClass.h: a copy between the kernels leaves the device idle for
// about 10 us, which was once for every launch sequence), each launch then reads its blocks from the offset that add returned
typedef struct BatchStageType6_ {
	std::vector<float> host;
	int blocks = 0;
	// Room for n blocks of the next launch sequence (zero), the offset of the first one in off. The pointer is only valid until the next add
	float* add(const int n, int& off) {
		off = blocks;
		host.resize(static_cast<size_t>(blocks + n) * type6BatchStride, 0.f);
		blocks += n;
		return &host[static_cast<size_t>(off) * type6BatchStride];
	}
} BatchStageType6;

// Whether the stacked arrays of nb views with the extents e stay below the size limit (type6BatchMaxElements), and the number of the
// groups along the third dimension of the launches below 65535
inline bool type6RunFits(const RunExtType6& e, const int nb) {
	const int nzMax = (std::max)(e.iz, e.nz);
	const int64_t slice = static_cast<int64_t>(roundUpType6(e.n0, type6AllocQuantum)) * roundUpType6(e.wy, type6AllocQuantum) * roundUpType6(nzMax, type6AllocQuantum);
	return static_cast<int64_t>(nb) * slice <= type6BatchMaxElements && static_cast<int64_t>(nb) * nzMax <= 65535;
}

// The parameter blocks (one for each view, type6BatchStride floats) of the forward projection of a run, see forwardRunType6
inline void type6FillRunFP(float* host, const std::vector<BatchItemType6>& run, const Weighting& w_vec) {
	const int nb = static_cast<int>(run.size());
	const bool resample = run[0].vv->grid.resample;
	for (int jj = 0; jj < nb; jj++) {
		const BatchItemType6& it = run[jj];
		float* p = &host[static_cast<size_t>(jj) * type6BatchStride];
		type6PutRotation(p, it.bwin, -w_vec.swivelAngles[it.u1], it.offX, it.offY);
		type6PutView(p, it);
		// The projection goes to the slot of the view in fProj, or, when it is resampled, to the slice of the array that the resampling reads
		p[type6P::OUTSLOT] = static_cast<float>(resample ? jj : it.kk);
	}
}

// Forward projection of volume ii for the views of a run (every one of them on the in-kernel weights, see type6RunFits and
// type6GlobalWeights): the rotation, the axial pass and the transaxial pass (and the resampling) are one launch each, and the
// projections are added to the projections of fProj. source is the volume, rotated through the image texSlot with images (it is
// uploaded when uploaded is not set). The attenuation is the one of the items (all of them use the same array or none). The parameter
// blocks of the run (type6FillRunFP) are in the device buffer of the batch, from the block blockOff on
inline void forwardRunType6(af::array& fProj, const std::vector<BatchItemType6>& run, const af::array& source, const int texSlot, char& uploaded, const Weighting& w_vec,
	const scalarStruct& inputScalars, ProjectorClass& proj, const int ii, const int blockOff) {
	const int nb = static_cast<int>(run.size());
	const int Q = type6AllocQuantum;
	const int nRowsD = static_cast<int>(inputScalars.nRowsD);
	const int nColsD = static_cast<int>(inputScalars.nColsD);
	const bool resample = run[0].vv->grid.resample;
	const bool useAtten = run[0].att != nullptr;
	RunExtType6 ext;
	for (int jj = 0; jj < nb; jj++)
		type6RunInclude(ext, run[jj]);
	if (!resample && (ext.ny != nRowsD || ext.nz != nColsD)) {
		mexPrint("The projection does not match the detector");
		return;
	}
	proj.prmOffBlocks = blockOff;
	// 1. Rotate the volume to the windows of the views
	if (!inputScalars.useBuffers && !uploaded) {
		if (rotUploadType6(source, run[0].win, inputScalars, proj, texSlot) != 0)
			return;
	}
	uploaded = 1;
	const int sXr = roundUpType6(ext.n0, Q), sYr = roundUpType6(ext.wy, Q), sZr = roundUpType6(ext.iz, Q);
	const int bsR = sXr * sYr * sZr;
	af::array rotated(af::dim4(sXr, sYr, sZr, nb), f32);
	proj.d_rhs = transferAF(rotated);
	if (inputScalars.useBuffers)
		proj.d_im = transferAF(source);
	int status = proj.rotateBatch(inputScalars, nb, ext.n0, ext.wy, ext.iz, sXr, sYr, sZr, texSlot);
	rotated.unlock();
	if (inputScalars.useBuffers)
		source.unlock();
	if (status != 0) {
		mexPrint("Failed to rotate the batch of views");
		return;
	}
	// 2. The axial pass (with the attenuation of the views)
	const int sXa = roundUpType6(ext.n0, Q), sYa = roundUpType6(ext.wy, Q), sZa = roundUpType6(ext.nz, Q);
	const int bsA = sXa * sYa * sZa;
	af::array axial(af::dim4(sXa, sYa, sZa, nb), f32);
	proj.d_im = transferAF(rotated);
	proj.d_rhs = transferAF(axial);
	int attSx = 0, attSy = 0;
	if (useAtten) {
		const af::array& attArr = *run[0].att->arr;
		proj.d_cdrfAtt = transferAF(attArr);
		attSx = static_cast<int>(run[0].att->arr->dims(0));
		attSy = static_cast<int>(run[0].att->arr->dims(1));
	}
	else
		proj.d_cdrfAtt = proj.d_im;
	status = proj.cdrfAxialBatch(inputScalars, nb, ext.n0, ext.wy, ext.nz, useAtten ? 1 : 0, run[0].prm, sXr, sYr, bsR, sXa, sYa, bsA, attSx, attSy);
	rotated.unlock();
	axial.unlock();
	if (useAtten)
		run[0].att->arr->unlock();
	if (status != 0) {
		mexPrint("Failed to apply the fused collimator response");
		return;
	}
	// 3. The transaxial pass, with the sum over the depth. Without resampling the projections are added to those of fProj
	if (!resample) {
		proj.d_im = transferAF(axial);
		proj.d_rhs = transferAF(fProj);
		status = proj.cdrfTransFPBatch(inputScalars, nb, ext.ny, ext.nz, sXa, sYa, bsA, nRowsD * nColsD, nRowsD, 1, run[0].prm);
		axial.unlock();
		fProj.unlock();
		if (status != 0)
			mexPrint("Failed to apply the fused collimator response");
		return;
	}
	// The projections on the lattice of the rotated grids are resampled onto the detector and added to fProj
	af::array summed(af::dim4(ext.ny, ext.nz, nb), f32);
	proj.d_im = transferAF(axial);
	proj.d_rhs = transferAF(summed);
	status = proj.cdrfTransFPBatch(inputScalars, nb, ext.ny, ext.nz, sXa, sYa, bsA, ext.ny * ext.nz, ext.ny, 0, run[0].prm);
	axial.unlock();
	summed.unlock();
	if (status == 0) {
		proj.d_im = transferAF(summed);
		proj.d_rhs = transferAF(fProj);
		status = proj.resampleToDetectorBatch(inputScalars, nb, nRowsD, nColsD, ext.ny, ext.ny * ext.nz, nRowsD * nColsD);
		summed.unlock();
		fProj.unlock();
	}
	if (status != 0)
		mexPrint("Failed to resample the projection onto the detector");
}

// The parameter blocks (one for each view, type6BatchStride floats) of the backprojection of a run, see backRunType6
inline void type6FillRunBP(float* host, const std::vector<BatchItemType6>& run, const Weighting& w_vec, const scalarStruct& inputScalars) {
	const int nb = static_cast<int>(run.size());
	const int nRowsD = static_cast<int>(inputScalars.nRowsD);
	const int nColsD = static_cast<int>(inputScalars.nColsD);
	const bool resample = run[0].vv->grid.resample;
	for (int jj = 0; jj < nb; jj++) {
		const BatchItemType6& it = run[jj];
		float* p = &host[static_cast<size_t>(jj) * type6BatchStride];
		type6PutRotation(p, it.bwin, w_vec.swivelAngles[it.u1], it.offX, it.offY);
		type6PutView(p, it);
		// The projection of the view: the measurements are resampled onto the lattice of the rotated grid inside the transaxial kernel, otherwise
		// the detector window is placed at the detector offset of the region (the measurements only cover the detector)
		if (resample) {
			p[type6P::PNY] = static_cast<float>(it.N1);
			p[type6P::PNZ] = static_cast<float>(it.N2);
			p[type6P::PY0] = 0.f;
			p[type6P::PZ0] = 0.f;
		}
		else {
			p[type6P::PNY] = static_cast<float>(nRowsD);
			p[type6P::PNZ] = static_cast<float>(nColsD);
			p[type6P::PY0] = it.vv->crop.cropped ? static_cast<float>(it.vv->crop.transOffset) : 0.f;
			p[type6P::PZ0] = it.vv->crop.cropped ? static_cast<float>(it.vv->crop.axialOffset) : 0.f;
		}
	}
}

// Backprojection of the projections fp of the views of a run of volume ii into output, the flat image of the volume: the spreading
// over the depth, the axial pass (with the attenuation) and the rotation back are one launch each, the last of them adds the views, in
// their order, to the image. The parameter blocks of the run (type6FillRunBP) are in the device buffer of the batch, from the block blockOff on
inline void backRunType6(af::array& output, const af::array& fp, const std::vector<BatchItemType6>& run, const Weighting& w_vec, const scalarStruct& inputScalars,
	ProjectorClass& proj, const int ii, const int blockOff) {
	const int nb = static_cast<int>(run.size());
	const int Q = type6AllocQuantum;
	const int nRowsD = static_cast<int>(inputScalars.nRowsD);
	const int nColsD = static_cast<int>(inputScalars.nColsD);
	const bool resample = run[0].vv->grid.resample;
	const bool useAtten = run[0].att != nullptr;
	const int nvx = static_cast<int>(inputScalars.Nx[ii]), nvy = static_cast<int>(inputScalars.Ny[ii]), nvz = static_cast<int>(inputScalars.Nz[ii]);
	if (static_cast<int64_t>(output.elements()) != static_cast<int64_t>(nvx) * nvy * nvz) {
		mexPrint("The rotation accumulates into an array of the wrong size");
		return;
	}
	RunExtType6 ext;
	for (int jj = 0; jj < nb; jj++)
		type6RunInclude(ext, run[jj]);
	proj.prmOffBlocks = blockOff;
	// The spread projections are on the axial rows of the input window of the axial pass, the regions on those of its output
	const int sX = roundUpType6(ext.n0, Q), sY = roundUpType6(ext.wy, Q), sZs = roundUpType6(ext.iz, Q), sZ = roundUpType6(ext.nz, Q);
	const int bsS = sX * sY * sZs;
	const int bs = sX * sY * sZ;
	// 1. Spread the projections over the depth and apply the transaxial response
	af::array spread(af::dim4(sX, sY, sZs, nb), f32);
	proj.d_im = transferAF(fp);
	proj.d_rhs = transferAF(spread);
	int status = proj.cdrfTransBPBatch(inputScalars, nb, ext.n0, ext.wy, ext.iz, sX, sY, bsS, run[0].prm, resample ? 1 : 0, nRowsD, nColsD);
	fp.unlock();
	spread.unlock();
	if (status != 0) {
		mexPrint("Failed to apply the fused collimator response");
		return;
	}
	// 2. The axial response and the attenuation
	af::array region(af::dim4(sX, sY, sZ, nb), f32);
	proj.d_im = transferAF(spread);
	proj.d_rhs = transferAF(region);
	int attSx = 0, attSy = 0;
	if (useAtten) {
		const af::array& attArr = *run[0].att->arr;
		proj.d_cdrfAtt = transferAF(attArr);
		attSx = static_cast<int>(run[0].att->arr->dims(0));
		attSy = static_cast<int>(run[0].att->arr->dims(1));
	}
	else
		proj.d_cdrfAtt = proj.d_im;
	status = proj.cdrfAxialBatch(inputScalars, nb, ext.n0, ext.wy, ext.nz, useAtten ? 2 : 0, run[0].prm, sX, sY, bsS, sX, sY, bs, attSx, attSy);
	spread.unlock();
	region.unlock();
	if (useAtten)
		run[0].att->arr->unlock();
	if (status != 0) {
		mexPrint("Failed to apply the fused collimator response");
		return;
	}
	// 3. Rotate the regions back into the image
	proj.d_im = transferAF(region);
	proj.d_rhs = transferAF(output);
	status = proj.rotateBackBatch(inputScalars, nb, nvx, nvy, nvz, sX, sY, bs);
	region.unlock();
	output.unlock();
	if (status != 0)
		mexPrint("Failed to rotate the batch of views back");
}

// The window of the rotation of the attenuation of group grp for one view (the views of vvb are those of the volumes), the union of those of the
// active members, see attenuationGroupType6. The rows and the axial rows of a member are those that its
// batched kernels read (type6AttWindowType6), the members whose view is empty are left out of them. The depth is not limited by them: the sum starts at the
// detector, so it is the same as without them. The function returns false when none of the members reads anything
typedef struct AttViewType6_ {
	int b = 0, first = -1;
	int rotX = 0, y0 = 0, y1 = 0, z0 = 0, z1 = 0;
	float offX = 0.f, offY = 0.f;
	RotWindow win;
} AttViewType6;

inline bool attViewWindowType6(AttViewType6& a, const AttGroupType6& grp, const std::vector<ViewVolumeType6>& vv, const int64_t u1, const Weighting& w_vec,
	const scalarStruct& inputScalars) {
	bool any = false;
	for (size_t mm = 0; mm < grp.members.size(); mm++) {
		const int ii = grp.members[mm];
		if (!vv[ii].active)
			continue;
		const int xoff = (grp.attNx - vv[ii].grid.Nx) / 2;
		const int rowOff = (grp.rows - vv[ii].grid.Ny) / 2;
		const int zShift = grp.zOff[ii] + grp.padLo - padZLoType6(inputScalars, ii);
		const int xB = xoff + vv[ii].dw.xEnd;
		if (a.first < 0) {
			a.first = ii;
			a.rotX = xB;
		}
		else
			a.rotX = (std::max)(a.rotX, xB);
		int wy0 = 0, wy1 = 0, wz0 = 0, wz1 = 0;
		if (!type6AttWindowType6(vv[ii], inputScalars, ii, wy0, wy1, wz0, wz1))
			continue;
		const int yA = vv[ii].crop.transStart + rowOff + wy0, yB = vv[ii].crop.transStart + rowOff + wy1;
		const int zA = vv[ii].crop.axialStart + zShift + wz0, zB = vv[ii].crop.axialStart + zShift + wz1;
		if (!any) {
			any = true;
			a.y0 = yA; a.y1 = yB; a.z0 = zA; a.z1 = zB;
		}
		else {
			a.y0 = (std::min)(a.y0, yA);
			a.y1 = (std::max)(a.y1, yB);
			a.z0 = (std::min)(a.z0, zA);
			a.z1 = (std::max)(a.z1, zB);
		}
	}
	if (a.first < 0 || !any)
		return false;
	RotWindow& win = a.win;
	win.srcNx = grp.srcNx; win.srcNyFull = grp.srcNy; win.srcNzFull = grp.padLo + grp.nz + grp.padHi;
	win.outNx = grp.attNx; win.outNyFull = grp.rows;
	win.outX = 0; win.outNxW = a.rotX;
	win.outY = a.y0; win.outNy = a.y1 - a.y0;
	win.outZ = a.z0; win.outNz = a.z1 - a.z0;
	win.srcY = 0; win.srcNy = grp.srcNy;
	win.srcZ = grp.padLo; win.srcNz = grp.nz;
	// The source is centered on the rotation axis, so it has no offset from it
	GridType6 attGrid = vv[a.first].grid;
	attGrid.axisX = 0.f;
	attGrid.axisY = 0.f;
	offsetsFPType6(w_vec, attGrid, u1, a.offX, a.offY);
	return true;
}

// What the launch of the attenuation of one group of a batch needs, made by attenuationGroupPlanType6: where its parameter blocks are in the stage
// of the batch (off, one for each of the na views that have an active member), the extents of the stack and the window of its first view
typedef struct AttLaunchType6_ {
	size_t gg = 0;
	int off = 0, na = 0;
	int mX = 0, mY = 0, mZ = 0;
	int sX = 0, sY = 0, sZ = 0;
	RotWindow win;
} AttLaunchType6;

// The attenuation of the volumes of attenuation group grp for the views of a batch (attenuationGroupType6 for every one of them): the
// attenuation source of the group is rotated to the windows of the views with one launch, into a stack of one slice for each view that has
// an active member, and accumulated along the depth (the sum without the voxel size, which is not exponentiated, see AttRefType6) with
// one af::accum of the stack, into attAcc. Every active volume of the group gets in attRefs[view][volume] the offset of its window in
// the slice of the view. attAcc has to stay where it is as long as attRefs is used.
// This is made in two steps: the planning (this function) only works on the host: it puts the parameter blocks of the rotation into the stage
// of the batch (see BatchStageType6), fills attRefs and adds the launch to pending; attenuationGroupLaunchType6 makes it. Returns false
// when no volume of the group is active in any of the views (nothing is planned)
inline bool attenuationGroupPlanType6(std::vector<AttLaunchType6>& pending, std::vector<std::vector<AttRefType6>>& attRefs, af::array& attAcc, const size_t gg, const AttGroupType6& grp,
	const std::vector<std::vector<ViewVolumeType6>>& vvb, const std::vector<int64_t>& u1s, const Weighting& w_vec, const scalarStruct& inputScalars, BatchStageType6& stage) {
	const int Q = type6AllocQuantum;
	// The windows of the views: the union of those of the active members, see attViewWindowType6
	std::vector<AttViewType6> av;
	for (size_t bb = 0; bb < vvb.size(); bb++) {
		AttViewType6 a;
		a.b = static_cast<int>(bb);
		if (attViewWindowType6(a, grp, vvb[bb], u1s[bb], w_vec, inputScalars))
			av.push_back(a);
	}
	if (av.empty())
		return false;
	const int na = static_cast<int>(av.size());
	int mX = 0, mY = 0, mZ = 0;
	AttLaunchType6 launch;
	float* host = stage.add(na, launch.off);
	for (int ss = 0; ss < na; ss++) {
		const AttViewType6& a = av[ss];
		mX = (std::max)(mX, a.win.outNxW);
		mY = (std::max)(mY, a.win.outNy);
		mZ = (std::max)(mZ, a.win.outNz);
		type6PutRotation(&host[static_cast<size_t>(ss) * type6BatchStride], a.win, -w_vec.swivelAngles[u1s[a.b]], a.offX, a.offY);
	}
	launch.gg = gg;
	launch.na = na;
	launch.mX = mX; launch.mY = mY; launch.mZ = mZ;
	launch.sX = roundUpType6(mX, Q); launch.sY = roundUpType6(mY, Q); launch.sZ = roundUpType6(mZ, Q);
	launch.win = av[0].win;
	pending.push_back(launch);
	for (int ss = 0; ss < na; ss++) {
		const AttViewType6& a = av[ss];
		for (size_t mm = 0; mm < grp.members.size(); mm++) {
			const int ii = grp.members[mm];
			const ViewVolumeType6& vv = vvb[a.b][ii];
			int wy0 = 0, wy1 = 0, wz0 = 0, wz1 = 0;
			if (!vv.active || !type6AttWindowType6(vv, inputScalars, ii, wy0, wy1, wz0, wz1))
				continue;
			const int xoff = (grp.attNx - vv.grid.Nx) / 2;
			const int rowOff = (grp.rows - vv.grid.Ny) / 2;
			const int zShift = grp.zOff[ii] + grp.padLo - padZLoType6(inputScalars, ii);
			AttRefType6& ref = attRefs[a.b][ii];
			ref.arr = &attAcc;
			ref.x0 = xoff + vv.dw.xOff;
			ref.y0 = vv.crop.transStart + rowOff - a.y0;
			ref.z0 = vv.crop.axialStart + zShift - a.z0;
			// The attenuation factors are exp(-dx * the accumulated attenuation), made by the kernel
			ref.dx = inputScalars.dx[a.first];
			ref.slice = ss;
			// The accumulated array has the dimensions of the stack
			ref.zStride = launch.sZ;
		}
	}
	return true;
}

// The launches of attenuationGroupPlanType6: the rotation of the attenuation source attSrc of the group (in the image texSlot, uploaded is set
// when it is already there) and the accumulation along the depth, into attAcc. baseOff is the block of the device buffer that the stage of the
// batch was uploaded to. The elements outside of the windows are not initialized, which does not matter, see attenuationGroupType6: the sum is
// along the first dimension, and each element only depends on the ones before it in its own line
inline bool attenuationGroupLaunchType6(af::array& attAcc, const AttLaunchType6& launch, const int baseOff, const af::array& attSrc, const scalarStruct& inputScalars,
	ProjectorClass& proj, const int texSlot, char& uploaded) {
	proj.prmOffBlocks = baseOff + launch.off;
	if (!inputScalars.useBuffers && !uploaded) {
		if (rotUploadType6(attSrc, launch.win, inputScalars, proj, texSlot) != 0)
			return false;
	}
	uploaded = 1;
	af::array attenuationImage(af::dim4(launch.sX, launch.sY, launch.sZ, launch.na), f32);
	proj.d_rhs = transferAF(attenuationImage);
	if (inputScalars.useBuffers)
		proj.d_im = transferAF(attSrc);
	const int status = proj.rotateBatch(inputScalars, launch.na, launch.mX, launch.mY, launch.mZ, launch.sX, launch.sY, launch.sZ, texSlot);
	attenuationImage.unlock();
	if (inputScalars.useBuffers)
		attSrc.unlock();
	if (status != 0) {
		mexPrint("Failed to rotate the attenuation of the batch of views");
		return false;
	}
	// Accumulate the attenuation along the depth
	attAcc = af::accum(attenuationImage, 0);
	return true;
}

// Everything that the attenuation of a batch of views (attenuationBatchPlanType6) is made from, besides the source and the groups (which are the key of
// the cache): for every view the projection and its angle, and for every group the windows of the active volumes and the offsets of the rotation
inline std::vector<double> attBatchSigType6(const std::vector<AttGroupType6>& groups, const std::vector<std::vector<ViewVolumeType6>>& vvb, const std::vector<int64_t>& u1s,
	const Weighting& w_vec, const scalarStruct& inputScalars) {
	std::vector<double> sig;
	for (size_t bb = 0; bb < vvb.size(); bb++) {
		sig.push_back(static_cast<double>(u1s[bb]));
		sig.push_back(w_vec.swivelAngles[u1s[bb]]);
		for (size_t gg = 0; gg < groups.size(); gg++) {
			const AttGroupType6& grp = groups[gg];
			bool first = true;
			for (size_t mm = 0; mm < grp.members.size(); mm++) {
				const int ii = grp.members[mm];
				const ViewVolumeType6& vv = vvb[bb][ii];
				sig.push_back(vv.active ? 1.0 : 0.0);
				if (!vv.active)
					continue;
				sig.push_back(vv.grid.Nx); sig.push_back(vv.grid.Ny);
				sig.push_back(vv.dw.xOff); sig.push_back(vv.dw.xEnd);
				sig.push_back(vv.crop.transStart); sig.push_back(vv.crop.transLength);
				sig.push_back(vv.crop.axialStart); sig.push_back(vv.crop.axialLength);
				sig.push_back(inputScalars.dx[ii]);
				int wy0 = 0, wy1 = 0, wz0 = 0, wz1 = 0;
				const bool reads = type6AttWindowType6(vv, inputScalars, ii, wy0, wy1, wz0, wz1);
				sig.push_back(reads ? 1.0 : 0.0);
				sig.push_back(wy0); sig.push_back(wy1); sig.push_back(wz0); sig.push_back(wz1);
				if (first) {
					// The offsets of the rotation, from the grid of the first active volume as attenuationGroupBatchType6 does
					GridType6 attGrid = vv.grid;
					attGrid.axisX = 0.f;
					attGrid.axisY = 0.f;
					float offX = 0.f, offY = 0.f;
					offsetsFPType6(w_vec, attGrid, u1s[bb], offX, offY);
					sig.push_back(offX); sig.push_back(offY);
					first = false;
				}
			}
			// The window of the rotation of the group (attViewWindowType6)
			AttViewType6 a;
			const bool made = attViewWindowType6(a, grp, vvb[bb], u1s[bb], w_vec, inputScalars);
			sig.push_back(made ? 1.0 : 0.0);
			if (made) {
				sig.push_back(a.first); sig.push_back(a.rotX);
				sig.push_back(a.y0); sig.push_back(a.y1); sig.push_back(a.z0); sig.push_back(a.z1);
				sig.push_back(a.offX); sig.push_back(a.offY);
			}
		}
	}
	return sig;
}

// The attenuation of a batch of views, in an AttBatchType6 (whose arrays and attRefs stay valid as long as it is held), planned: nothing is
// launched, the parameter blocks of the rotations are in the stage of the batch and the launches in pending (attenuationBatchLaunchType6 makes
// them after the stage has been uploaded). With the cache (cache is not null): the forward projection (store) keeps the batch for the
// backprojection of its subset, as long as the kept batches stay within budget bytes; the backprojection takes the batch that the forward
// projection kept when it is the same one (the same views, the same source and windows), and removes it from the cache (there is nothing to
// launch then), otherwise makes it as without the cache. The kept arrays are exactly those that the backprojection would make
inline std::shared_ptr<AttBatchType6> attenuationBatchPlanType6(const bool store, AttCacheType6* cache, const bool useAtten, const std::vector<AttGroupType6>& groups,
	const std::vector<std::vector<ViewVolumeType6>>& vvb, const std::vector<int64_t>& u1s, const int64_t uu, const int64_t length, const int64_t k0, const int nb, const int nV,
	const Weighting& w_vec, const scalarStruct& inputScalars, BatchStageType6& stage, std::vector<AttLaunchType6>& pending) {
	std::shared_ptr<AttBatchType6> batch = std::make_shared<AttBatchType6>();
	batch->refs.assign(static_cast<size_t>(nb), std::vector<AttRefType6>(static_cast<size_t>(nV)));
	batch->acc.assign(groups.size(), af::array());
	if (!useAtten)
		return batch;
	if (cache != nullptr) {
		batch->uu = uu; batch->length = length; batch->k0 = k0; batch->nb = nb;
		batch->sig = attBatchSigType6(groups, vvb, u1s, w_vec, inputScalars);
		// The attenuation of a batch of views is made from the source of the map of the frame (see attenuationFrameType6)
		batch->sig.push_back(static_cast<double>(cache->frame));
		if (!store) {
			for (size_t cc = 0; cc < cache->batches.size(); cc++) {
				const std::shared_ptr<AttBatchType6> kept = cache->batches[cc];
				if (kept->uu == uu && kept->length == length && kept->k0 == k0 && kept->nb == nb && kept->acc.size() == groups.size() && kept->sig == batch->sig) {
					cache->batchBytes -= kept->bytes;
					cache->batches.erase(cache->batches.begin() + static_cast<std::ptrdiff_t>(cc));
					return kept;
				}
			}
		}
	}
	for (size_t gg = 0; gg < groups.size(); gg++)
		attenuationGroupPlanType6(pending, batch->refs, batch->acc[gg], gg, groups[gg], vvb, u1s, w_vec, inputScalars, stage);
	return batch;
}

// The launches planned by attenuationBatchPlanType6 for the batch, one group at a time (baseOff: where the stage of the batch is in the device
// buffer), and, with the cache, the keeping of the batch (the forward projection: store) when it fits within the budget. Nothing is done for
// a batch that the backprojection took from the cache (pending is empty and store is not set)
inline void attenuationBatchLaunchType6(const std::shared_ptr<AttBatchType6>& batch, const std::vector<AttLaunchType6>& pending, const bool store, AttCacheType6* cache, const int64_t budget,
	const bool useAtten, const int baseOff, const std::vector<af::array>& attenSrc, std::vector<char>& attUploaded, const scalarStruct& inputScalars, ProjectorClass& proj) {
	if (!useAtten)
		return;
	for (size_t pp = 0; pp < pending.size(); pp++) {
		const size_t gg = pending[pp].gg;
		attenuationGroupLaunchType6(batch->acc[gg], pending[pp], baseOff, attenSrc[gg], inputScalars, proj, rotSlotAttenuationType6 + static_cast<int>(gg), attUploaded[gg]);
	}
	if (cache != nullptr && store) {
		int64_t bytes = 0;
		for (size_t gg = 0; gg < batch->acc.size(); gg++)
			bytes += static_cast<int64_t>(batch->acc[gg].elements()) * static_cast<int64_t>(sizeof(float));
		if (cache->batchBytes + bytes <= budget) {
			batch->bytes = bytes;
			cache->batches.push_back(batch);
			cache->batchBytes += bytes;
		}
	}
}

// One operation of a batch of views on one volume, in the order that the batch makes them in: a run of views (off is the block of its parameters in
// the stage of the batch) or, when off is negative, a view on the global weights that is processed alone (the item)
typedef struct BatchOpType6_ {
	int ii = 0;
	int off = -1;
	std::vector<BatchItemType6> run;
	BatchItemType6 item;
} BatchOpType6;

// SPECT forward projection (projector type 6) of every volume, with several views per kernel launch (see above and type6BatchSize). The
// arguments and the result are those of forwardProjectionType6All
inline void forwardProjectionType6AllBatched(af::array& fProj, const Weighting& w_vec, AF_im_vectors& vec, const scalarStruct& inputScalars,
	const int64_t length, const int64_t uu, ProjectorClass& proj, const float* atten, const int timestep) {
	if (DEBUG || inputScalars.verbose >= 3) {
		af::sync();
		proj.tStartLocal = std::chrono::steady_clock::now();
	}
	if (DEBUG || inputScalars.verbose >= 3)
		mexPrint("Starting SPECT forward projection of all volumes, several views per launch");
	const int nV = static_cast<int>(inputScalars.nMultiVolumes) + 1;
	const int bMax = type6BatchSizeEffective(inputScalars);
	const bool useAtten = (inputScalars.attenuation_correction && inputScalars.CTAttenuation && (atten != nullptr));
	std::vector<af::array> apuArr(nV);
	for (int ii = 0; ii < nV; ii++)
		apuArr[ii] = ownedType6(af::moddims(vec.im_os[timestep][ii], inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii]));
	const std::vector<AttGroupType6> groups = attGroupsType6(inputScalars);
	AttSourcesType6 attSources;
	attenuationSourcesCachedType6(attSources, proj, atten, inputScalars, groups, useAtten, static_cast<int>(timestep));
	const std::vector<af::array>& attenSrc = *attSources.src;
	std::vector<char>& attUploaded = *attSources.uploaded;
	std::vector<char> emUploaded(nV, 0);
	// The attenuation of the batches of views of this subset is kept for the backprojection of the subset (see attenuationBatchSharedType6);
	// what the previous subset kept is dropped
	AttCacheType6* attCache = useAtten ? &attCacheType6(proj) : nullptr;
	int64_t attBudget = 0;
	if (attCache != nullptr) {
		attCache->batches.clear();
		attCache->batchBytes = 0;
		attBudget = attBudgetType6(proj, inputScalars);
	}
	// The parameter blocks of all the launches of this call are in one device buffer (see beginBatchParams): at most one block for each view and
	// volume, and one for each view and attenuation group
	if (proj.beginBatchParams(static_cast<int>(length) * (nV + static_cast<int>(groups.size()))) != 0)
		return;
	for (int64_t k0 = 0; k0 < length; k0 += bMax) {
		const int nb = static_cast<int>((std::min)(static_cast<int64_t>(bMax), length - k0));
		std::vector<std::vector<ViewVolumeType6>> vvb(nb, std::vector<ViewVolumeType6>(nV));
		std::vector<int64_t> u1s(nb);
		for (int jj = 0; jj < nb; jj++) {
			u1s[jj] = uu + geomViewStartType6(w_vec, static_cast<int64_t>(timestep)) + k0 + jj;
			for (int ii = 0; ii < nV; ii++)
				vvb[jj][ii] = viewVolumeType6(w_vec, inputScalars, u1s[jj], ii);
		}
		// The batch is planned first (the parameter blocks of the attenuation and of every run go into one stage, which is uploaded with one copy),
		// then made
		BatchStageType6 stage;
		std::vector<AttLaunchType6> pendingAtt;
		const std::shared_ptr<AttBatchType6> attBatch = attenuationBatchPlanType6(true, attCache, useAtten, groups, vvb, u1s, uu, length,
			k0, nb, nV, w_vec, inputScalars, stage, pendingAtt);
		std::vector<std::vector<AttRefType6>>& attRefs = attBatch->refs;
		std::vector<BatchOpType6> ops;
		for (int ii = 0; ii < nV; ii++) {
			std::vector<BatchItemType6> run;
			RunExtType6 ext;
			const auto flush = [&]() {
				if (!run.empty()) {
					BatchOpType6 op;
					op.ii = ii;
					type6FillRunFP(stage.add(static_cast<int>(run.size()), op.off), run, w_vec);
					op.run = run;
					ops.push_back(op);
				}
				run.clear();
				ext = RunExtType6();
			};
			for (int jj = 0; jj < nb; jj++) {
				ViewVolumeType6& vv = vvb[jj][ii];
				if (!vv.active)
					continue;
				BatchItemType6 it;
				type6ItemFP(it, vv, useAtten ? &attRefs[jj][ii] : nullptr, static_cast<int>(k0) + jj, u1s[jj], w_vec, inputScalars, ii);
				if (it.empty)
					continue;
				if (it.prm.globalW) {
					// The global weights: this view is processed alone, as without batching
					BatchOpType6 op;
					op.ii = ii;
					op.item = it;
					ops.push_back(op);
					continue;
				}
				RunExtType6 t = ext;
				type6RunInclude(t, it);
				if (!run.empty() && (static_cast<int>(run.size()) >= bMax || !type6RunFits(t, static_cast<int>(run.size()) + 1))) {
					flush();
					t = RunExtType6();
					type6RunInclude(t, it);
				}
				run.push_back(it);
				ext = t;
			}
			flush();
		}
		int baseOff = 0;
		if (stage.blocks > 0) {
			if (proj.uploadBatchParams(stage.host.data(), stage.blocks) != 0)
				return;
			baseOff = proj.prmOffBlocks;
		}
		attenuationBatchLaunchType6(attBatch, pendingAtt, true, attCache, attBudget, useAtten, baseOff, attenSrc, attUploaded, inputScalars, proj);
		for (size_t oo = 0; oo < ops.size(); oo++) {
			const BatchOpType6& op = ops[oo];
			const int ii = op.ii;
			if (op.off >= 0) {
				forwardRunType6(fProj, op.run, apuArr[ii], rotSlotEmissionType6 + ii, emUploaded[ii], w_vec, inputScalars, proj, ii, baseOff + op.off);
				continue;
			}
			const BatchItemType6& it = op.item;
			ViewVolumeType6& vv = *it.vv;
			RotWindow rotWin = it.win;
			rotWin.allocQ = type6AllocQuantum;
			af::array kuvaRot = rotateHelperType6(apuArr[ii], inputScalars, proj, -w_vec.swivelAngles[it.u1], it.offX, it.offY, ii, rotSlotEmissionType6 + ii, !emUploaded[ii], rotWin);
			emUploaded[ii] = 1;
			prepareFusedType6(vv.PSF, vv.grid, inputScalars, proj, ii, vv.dw, vv.crop);
			forwardViewFusedType6(fProj, it.kk, kuvaRot, it.att, vv.grid, vv.PSF, vv.crop, vv.map, inputScalars, proj, ii, vv.dw, type6AllocQuantum);
		}
		if (DEBUG || inputScalars.verbose > 2)
			mexPrint("Projector 6 FP of all volumes, batch of views complete");
	}
	if (DEBUG || inputScalars.verbose >= 3) {
		proj.tEndLocal = std::chrono::steady_clock::now();
		const std::chrono::duration<double> tDiff = proj.tEndLocal - proj.tStartLocal;
		mexPrintBase("SPECT forward projection of all volumes completed in %f seconds\n", tDiff);
		mexEval();
	}
}

// SPECT backprojection (projector type 6) of every volume, with several views per kernel launch (see above and type6BatchSize). The
// arguments and the result are those of backprojectionType6All
inline void backprojectionType6AllBatched(af::array& fProj, const Weighting& w_vec, AF_im_vectors& vec,
	const scalarStruct& inputScalars, const int64_t length, const int64_t uu, ProjectorClass& proj, uint32_t timestep, const uint32_t osa_iter,
	const uint32_t iter, const uint8_t compute_norm_matrix, const uint32_t iter0, const float* atten, const uint32_t sensIdx, const bool allowSensCompute) {
	if (DEBUG || inputScalars.verbose >= 3)
		mexPrint("Starting SPECT backprojection of all volumes, several views per launch");
	const int nV = static_cast<int>(inputScalars.nMultiVolumes) + 1;
	const int bMax = type6BatchSizeEffective(inputScalars);
	const bool compSens = allowSensCompute && ((iter == iter0 && compute_norm_matrix == 2) || compute_norm_matrix == 1);
	af::array sensProj;
	// The adjoint of the forward projection of the model y ~ n*s*g*(A x): A^T (w .* data), and the sensitivity A^T w (see type6Weights)
	const af::array wgt = type6Weights(inputScalars, w_vec, uu, length, static_cast<int>(timestep));
	const af::array fIn = wgt.isempty() ? fProj : af::array(af::flat(fProj) * af::flat(wgt));
	if (compSens) {
		if (DEBUG || inputScalars.verbose >= 3)
			mexPrint("Computing sensitivity image");
		sensProj = wgt.isempty() ? af::constant(1.f, inputScalars.nRowsD, inputScalars.nColsD, length) : wgt;
	}
	const af::array fp = ownedType6(af::moddims(fIn, inputScalars.nRowsD, inputScalars.nColsD, length));
	af::array fpSens;
	if (compSens)
		fpSens = ownedType6(af::moddims(sensProj, inputScalars.nRowsD, inputScalars.nColsD, length));
	std::vector<af::array> output(nV), outputSens(nV);
	for (int ii = 0; ii < nV; ii++) {
		output[ii] = af::constant(0.f, inputScalars.Nx[ii] * inputScalars.Ny[ii] * inputScalars.Nz[ii], 1);
		if (compSens)
			outputSens[ii] = af::constant(0.f, inputScalars.Nx[ii] * inputScalars.Ny[ii] * inputScalars.Nz[ii], 1);
	}
	const bool useAtten = (inputScalars.attenuation_correction && inputScalars.CTAttenuation && (atten != nullptr));
	const std::vector<AttGroupType6> groups = attGroupsType6(inputScalars);
	AttSourcesType6 attSources;
	attenuationSourcesCachedType6(attSources, proj, atten, inputScalars, groups, useAtten, static_cast<int>(timestep));
	const std::vector<af::array>& attenSrc = *attSources.src;
	std::vector<char>& attUploaded = *attSources.uploaded;
	// The attenuation that the forward projection of this subset kept is used when it is the same (see attenuationBatchSharedType6)
	AttCacheType6* attCache = useAtten ? &attCacheType6(proj) : nullptr;
	const AttRefType6 noAtten;
	// The parameter blocks of all the launches of this call are in one device buffer (see beginBatchParams): at most one block for each view and
	// volume (twice with the sensitivity) and one for each view and attenuation group
	if (proj.beginBatchParams(static_cast<int>(length) * ((compSens ? 2 : 1) * nV + static_cast<int>(groups.size()))) != 0)
		return;
	for (int64_t k0 = 0; k0 < length; k0 += bMax) {
		const int nb = static_cast<int>((std::min)(static_cast<int64_t>(bMax), length - k0));
		std::vector<std::vector<ViewVolumeType6>> vvb(nb, std::vector<ViewVolumeType6>(nV));
		std::vector<int64_t> u1s(nb);
		for (int jj = 0; jj < nb; jj++) {
			u1s[jj] = uu + geomViewStartType6(w_vec, static_cast<int64_t>(timestep)) + k0 + jj;
			for (int ii = 0; ii < nV; ii++)
				vvb[jj][ii] = viewVolumeType6(w_vec, inputScalars, u1s[jj], ii);
		}
		// The same transformation as in the forward projection, accumulated along the projection direction (the kernel exponentiates). The batch is
		// planned first (one stage of parameter blocks, one copy), then made, as in the forward projection
		BatchStageType6 stage;
		std::vector<AttLaunchType6> pendingAtt;
		const std::shared_ptr<AttBatchType6> attBatch = attenuationBatchPlanType6(false, attCache, useAtten, groups, vvb, u1s, uu, length,
			k0, nb, nV, w_vec, inputScalars, stage, pendingAtt);
		std::vector<std::vector<AttRefType6>>& attRefs = attBatch->refs;
		std::vector<BatchOpType6> ops;
		for (int ii = 0; ii < nV; ii++) {
			std::vector<BatchItemType6> run;
			RunExtType6 ext;
			const auto flush = [&]() {
				if (!run.empty()) {
					BatchOpType6 op;
					op.ii = ii;
					type6FillRunBP(stage.add(static_cast<int>(run.size()), op.off), run, w_vec, inputScalars);
					op.run = run;
					ops.push_back(op);
				}
				run.clear();
				ext = RunExtType6();
			};
			for (int jj = 0; jj < nb; jj++) {
				ViewVolumeType6& vv = vvb[jj][ii];
				if (!vv.active)
					continue;
				BatchItemType6 it;
				type6ItemBP(it, vv, useAtten ? &attRefs[jj][ii] : nullptr, static_cast<int>(k0) + jj, u1s[jj], w_vec, inputScalars, ii);
				if (it.empty)
					continue;
				if (it.prm.globalW) {
					// The global weights: this view is processed alone, as without batching. The views that precede it have to be in the image first,
					// which keeps the order of the sums
					flush();
					BatchOpType6 op;
					op.ii = ii;
					op.item = it;
					ops.push_back(op);
					continue;
				}
				RunExtType6 t = ext;
				type6RunInclude(t, it);
				if (!run.empty() && (static_cast<int>(run.size()) >= bMax || !type6RunFits(t, static_cast<int>(run.size()) + 1))) {
					flush();
					t = RunExtType6();
					type6RunInclude(t, it);
				}
				run.push_back(it);
				ext = t;
			}
			flush();
		}
		int baseOff = 0;
		if (stage.blocks > 0) {
			if (proj.uploadBatchParams(stage.host.data(), stage.blocks) != 0)
				return;
			baseOff = proj.prmOffBlocks;
		}
		attenuationBatchLaunchType6(attBatch, pendingAtt, false, attCache, 0, useAtten, baseOff, attenSrc, attUploaded, inputScalars, proj);
		for (size_t oo = 0; oo < ops.size(); oo++) {
			const BatchOpType6& op = ops[oo];
			const int ii = op.ii;
			if (op.off >= 0) {
				// The sensitivity uses the parameters of the run as they are
				backRunType6(output[ii], fp, op.run, w_vec, inputScalars, proj, ii, baseOff + op.off);
				if (compSens)
					backRunType6(outputSens[ii], fpSens, op.run, w_vec, inputScalars, proj, ii, baseOff + op.off);
				continue;
			}
			const BatchItemType6& it = op.item;
			ViewVolumeType6& vv = *it.vv;
			prepareFusedType6(vv.PSF, vv.grid, inputScalars, proj, ii, vv.dw, vv.crop);
			const AttRefType6& attenuation = useAtten ? *it.att : noAtten;
			backProjectionType6Angle(output[ii], fp, it.kk, vv.grid, vv.PSF, vv.crop, vv.map, attenuation, useAtten, w_vec, inputScalars, proj, it.u1, ii, vv.dw,
				type6AllocQuantum);
			if (compSens)
				backProjectionType6Angle(outputSens[ii], fpSens, it.kk, vv.grid, vv.PSF, vv.crop, vv.map, attenuation, useAtten, w_vec, inputScalars, proj, it.u1, ii, vv.dw,
					type6AllocQuantum);
		}
	}
	for (int ii = 0; ii < nV; ii++) {
		// The mask of every volume (see prepareMaskBPType6)
		if (w_vec.maskBP != nullptr) {
			if (!applyMaskBP(output[ii], compSens ? &outputSens[ii] : nullptr, inputScalars, w_vec, ii, static_cast<int>(timestep)))
				return;
		}
		vec.rhs_os[timestep][ii] = output[ii];
		vec.rhs_os[timestep][ii](vec.rhs_os[timestep][ii] < inputScalars.epps && vec.rhs_os[timestep][ii] >= 0.f) = inputScalars.epps;
		if (compSens) {
			if (compute_norm_matrix == 2) {
				vec.Summ[sensIdx][ii][osa_iter] = outputSens[ii];
				vec.Summ[sensIdx][ii][osa_iter](vec.Summ[sensIdx][ii][osa_iter] < inputScalars.epps) = inputScalars.epps;
			}
			else {
				vec.Summ[sensIdx][ii][0] = outputSens[ii];
				vec.Summ[sensIdx][ii][0](vec.Summ[sensIdx][ii][0] < inputScalars.epps) = inputScalars.epps;
			}
			if (DEBUG || inputScalars.verbose >= 3)
				mexPrint("Sensitivity image computed");
		}
		if (w_vec.maskBP != nullptr) {
			if (!resetMaskedBP(vec.rhs_os[timestep][ii], compSens ? &vec.Summ[sensIdx][ii][compute_norm_matrix == 2 ? osa_iter : 0] : nullptr, inputScalars, w_vec, ii, static_cast<int>(timestep)))
				return;
		}
	}
	if (DEBUG || inputScalars.verbose >= 3)
		mexPrint("SPECT backprojection of all volumes complete");
}
#endif


// 1D filtering
inline int filtering(const af::array& filter, af::array& input, ProjectorClass& proj, const dim_t dimmi) {
	int status = 0;
	af::array temp = af::fft(input, dimmi);
	temp.eval();
	if (DEBUG) {
		mexPrintBase("temp = %f\n", af::sum<float>(temp));
		mexEval();
	}
#ifndef CPU
	status = elementWiseAF(filter, temp, true, proj);
	if (status != 0)
		return -1;
#else
	temp *= af::tile(filter, 1, temp.dims(1), temp.dims(2));
#endif
	if (DEBUG) {
		mexPrintBase("temp = %f\n", af::sum<float>(temp));
		mexEval();
	}
	af::ifftInPlace(temp);
	temp.eval();
	input = af::flat(af::real(temp(af::seq(0, input.dims(0) - 1), af::span, af::span)));
	return 0;
}

// 2D filtering
inline int filtering2D(const af::array& filter, af::array& input, ProjectorClass& proj, const dim_t dimmi) {
	int status = 0;
	const dim_t inDim0 = input.dims(0);
	if (DEBUG) {
		mexPrintBase("filter.dims(0) = %d\n", filter.dims(0));
		mexPrintBase("filter.dims(1) = %d\n", filter.dims(1));
		mexPrintBase("input.dims(0) = %d\n", input.dims(0));
		mexPrintBase("input.dims(1) = %d\n", input.dims(1));
		mexPrintBase("input.dims(2) = %d\n", input.dims(2));
		mexPrintBase("dimmi = %d\n", dimmi);
		mexEval();
	}
	af::array temp = af::fft2(input, dimmi, dimmi);
#ifndef CPU
	status = elementWiseAF(filter, temp, true, proj, true);
	if (status != 0)
		return -1;
#else
	temp *= filter;
#endif
	ifft2InPlace(temp);
	input = af::flat(af::real(temp(af::seq(0, input.dims(0) - 1), af::seq(0, input.dims(0) - 1), af::span)));
	input.eval();
	return 0;
}

// 1D filtering with division
inline int filteringInverse(const af::array& filter, af::array& input, ProjectorClass& proj, const dim_t dimmi) {
	int status = 0;
	af::array temp = af::fft(input, dimmi);
	temp.eval();
#ifndef CPU
	status = elementWiseAF(filter, temp, false, proj);
	if (status != 0)
		return -1;
#else
	temp /= af::tile(filter, 1, temp.dims(1), temp.dims(2));
	temp.eval();
#endif
	af::ifftInPlace(temp);
	input = af::flat(af::real(temp(af::seq(0, input.dims(0) - 1), af::span, af::span)));
	return 0;
}

// Image-based preconditioner 2
inline af::array precondIm2(const af::array& im, const af::array& D) {
	return im / D;
}

// Image-based preconditioner 3
inline af::array precondIm3(const af::array& im, const af::array& D, const af::array& ref) {
	return (af::max)(im, (af::max)(VAL, ref)) / D;
}

// Compute the gradient of an array
inline void computeGradient(const af::array& im, const scalarStruct& inputScalars, af::array& f, af::array& g, af::array& h, const int type = 0) {
	if (DEBUG) {
		mexPrintBase("im.dims(0) = %d\n", im.dims(0));
		mexPrintBase("im.dims(1) = %d\n", im.dims(1));
		mexPrintBase("im.dims(2) = %d\n", im.dims(2));
		mexEval();
	}
	// 1st order forward differences
	if (type == 0) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting forward difference gradient");
		f(af::seq(0, af::end - 1LL), af::span, af::span) = af::diff1(im);
		f(af::end, af::span, af::span) = -1.f * im(af::end, af::span, af::span);
		g(af::span, af::seq(0, af::end - 1LL), af::span) = af::diff1(im, 1);
		g(af::span, af::end, af::span) = -1.f * im(af::span, af::end, af::span);
		h(af::span, af::span, af::seq(0, af::end - 1LL)) = af::diff1(im, 2);
		h(af::span, af::span, af::end) = -1.f * im(af::span, af::span, af::end);
	}
	// 1st order backward differences
	else if (type == 1) {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting backward difference gradient");
		f(af::seq(1, af::end), af::span, af::span) = -af::diff1(im);
		f(0, af::span, af::span) = -1.f * im(0, af::span, af::span);
		g(af::span, af::seq(1, af::end), af::span) = -af::diff1(im, 1);
		g(af::span, 0, af::span) = -1.f * im(af::span, 0, af::span);
		h(af::span, af::span, af::seq(1, af::end)) = -af::diff1(im, 2);
		h(af::span, af::span, 0) = -1.f * im(af::span, af::span, 0);
	}
	// 1st order central differences
	else {
		if (inputScalars.verbose >= 3)
			mexPrint("Starting central difference gradient");
		f = (af::shift(im, -1) - af::shift(im, 1)) * .5f;
		f(0, af::span, af::span) = im(1, af::span, af::span) - im(0, af::span, af::span);
		f(af::end, af::span, af::span) = im(af::end, af::span, af::span) - im(af::end - 1, af::span, af::span);
		g = (af::shift(im, 0, -1) - af::shift(im, 0, 1)) * .5f;
		g(af::span, 0, af::span) = im(af::span, 1, af::span) - im(af::span, 0, af::span);
		g(af::span, af::end, af::span) = im(af::span, af::end, af::span) - im(af::span, af::end - 1, af::span);
		h = (af::shift(im, 0, 0, -1) - af::shift(im, 0, 0, 1)) * .5f;
		h(af::span, af::span, 0) = im(af::span, af::span, 1) - im(af::span, af::span, 0);
		h(af::span, af::span, af::end) = im(af::span, af::span, af::end) - im(af::span, af::span, af::end - 1);
	}
	f = af::flat(f);
	g = af::flat(g);
	h = af::flat(h);
	f.eval();
	g.eval();
	h.eval();
	if (inputScalars.verbose >= 3)
		mexPrint("Gradient computed");
}

// Second order gradient
inline void computeSecondOrderGradient(af::array& f, af::array& g, af::array& h, const int type = 0) {
	// Second order derivatives
	if (type == 0) {
		f = af::shift(f, 1);
		g = af::shift(g, 0, 1);
		h = af::shift(h, 0, 0, 1);
	}
	else if (type == 1) {
		f = af::shift(f, -1);
		g = af::shift(g, 0, -1);
		h = af::shift(h, 0, 0, -1);
	}
	else {
		f = (af::shift(f, -1) - af::shift(f, 1)) * .5f;
		g = (af::shift(g, 0, -1) - af::shift(g, 0, 1)) * .5f;
		h = (af::shift(h, 0, 0, -1) - af::shift(h, 0, 0, 1)) * .5f;
	}
}

// Image-based preconditioner 4
inline void gradientPreconditioner(const scalarStruct& inputScalars, Weighting& w_vec, const af::array& input, const int ii = 0) {
	af::array f = af::constant(0.f, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii]);
	af::array g = af::constant(0.f, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii]);
	af::array h = af::constant(0.f, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii]);
	computeGradient(input, inputScalars, f, g, h, w_vec.derivType);
	if (DEBUG) {
		mexPrintBase("g.dims(0) = %d\n", g.dims(0));
		mexEval();
	}
	f = (af::max)(VAL, af::sqrt(f * f + g * g + h * h) / af::mean<float>(af::flat(input)));
	f = af::mean<float>(af::flat(f)) / f;
	f.eval();
	if (DEBUG) {
		mexPrintBase("f.dims(0) = %d\n", f.dims(0));
		mexEval();
	}
	w_vec.gradF[ii] = (af::min)(w_vec.gradV2, (af::max)(w_vec.gradV1, f));
}

// Compute the image-based preconditioning
inline int applyImagePreconditioning(Weighting& w_vec, const scalarStruct& inputScalars, af::array& input, const af::array& im, ProjectorClass& proj, uint32_t timestep, const int kk = 0, const int ii = 0) {
	int status = 0;
	if (w_vec.precondTypeIm[4] && kk >= w_vec.gradInitIter) {
		if (inputScalars.verbose >= 3)
			mexPrint("Applying gradient-based preconditioner, type 4");
		if (kk <= w_vec.gradFinalIter)
			gradientPreconditioner(inputScalars, w_vec, af::moddims(im, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii]), ii);
		input *= w_vec.gradF[ii];
	}
	if (w_vec.precondTypeIm[3]) {
		if (inputScalars.verbose >= 3)
			mexPrint("Applying momentum-like preconditioner, type 3");
		input *= w_vec.alphaPrecond[kk];
	}
	if (w_vec.precondTypeIm[0] || w_vec.precondTypeIm[1] || w_vec.precondTypeIm[2]) {
		if (w_vec.precondTypeIm[0]) {
			if (inputScalars.verbose >= 3)
				mexPrint("Applying diagonal normalization preconditioner , type 0");
			input /= w_vec.D[0][ii];
		}
		else if (w_vec.precondTypeIm[1]) {
			if (inputScalars.verbose >= 3)
				mexPrint("Applying EM preconditioner, type 1");
			input *= precondIm2(im, w_vec.D[0][ii]);
		}
		else if (w_vec.precondTypeIm[2]) {
			if (inputScalars.verbose >= 3)
				mexPrint("Applying IEM preconditioner, type 2");
			input *= precondIm3(im, w_vec.D[0][ii], w_vec.preRef[ii]);
		}
	}
	if (w_vec.precondTypeIm[6]) {
		if (inputScalars.verbose >= 3)
			mexPrint("Applying curvature preconditioner , type 6");
		input *= w_vec.dP[timestep][ii];
	}
	if (w_vec.precondTypeIm[5] && kk <= w_vec.filterIter) {
		if (inputScalars.verbose >= 3)
			mexPrint("Applying filtering-based preconditioner, type 5");
		input = af::moddims(input, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii]);
		status = filtering2D(w_vec.filterIm, input, proj, inputScalars.Nf);
		if (status != 0)
			return -1;
	}
	input.eval();
	if (inputScalars.verbose >= 3 && (w_vec.precondTypeIm[0] || w_vec.precondTypeIm[1] || w_vec.precondTypeIm[2] || w_vec.precondTypeIm[3] ||
		(w_vec.precondTypeIm[4] && kk >= w_vec.gradInitIter) || w_vec.precondTypeIm[5] || w_vec.precondTypeIm[6]))
		mexPrint("Image-based preconditioning applied");
	return 0;
}

// Compute measurement-based preconditioning
inline int applyMeasPreconditioning(const Weighting& w_vec, const scalarStruct& inputScalars, af::array& input, ProjectorClass& proj, const uint32_t timestep, const uint32_t subIter = 0) {
	int status = 0;
	if (w_vec.precondTypeMeas[0] || w_vec.precondTypeMeas[1]) {
		if (inputScalars.verbose >= 3)
			mexPrint("Applying measurement-based preconditioning");
		if (w_vec.precondTypeMeas[1]) {
			if (inputScalars.verbose >= 3)
				mexPrint("Applying filtering-based preconditioner, type 1");
			if (DEBUG) {
				mexPrintBase("input.elements() = %d\n", input.elements());
				mexPrintBase("w_vec.filter.elements() = %d\n", w_vec.filter.elements());
				mexPrintBase("w_vec.filter = %f\n", af::sum<float>(w_vec.filter));
				mexPrintBase("inputScalars.nRowsD = %d\n", inputScalars.nRowsD);
				mexPrintBase("inputScalars.nColsD = %d\n", inputScalars.nColsD);
				mexPrintBase("input.elements() / (inputScalars.nRowsD * inputScalars.nColsD) = %d\n", input.elements() / (inputScalars.nRowsD * inputScalars.nColsD));
				mexEval();
			}
			if (inputScalars.subsetsUsed > 1 && (inputScalars.subsetType == 5 || inputScalars.subsetType == 4)) {
				if (inputScalars.subsetType == 4)
					input = af::moddims(input, inputScalars.nRowsD, input.elements() / inputScalars.nRowsD);
				else
					input = af::moddims(input, inputScalars.nColsD, input.elements() / inputScalars.nColsD);
			}
			else
				input = af::moddims(input, inputScalars.nRowsD, inputScalars.nColsD, input.elements() / (inputScalars.nRowsD * inputScalars.nColsD));
			input.eval();
			if (input.elements() > 207762000) {
				const dim_t h = input.dims(2) / 2;
				for (int uu = 0; uu < 2; uu++) {
					const dim_t first = uu * h, last = (uu == 0) ? h - 1 : input.dims(2) - 1;
					af::array input1 = input(af::span, af::span, af::seq(first, last)).copy();
					status = filtering(w_vec.filter, input1, proj, inputScalars.Nf);
					input1 = af::moddims(input1, inputScalars.nRowsD, inputScalars.nColsD, input1.elements() / (inputScalars.nRowsD * inputScalars.nColsD));
					input(af::span, af::span, af::seq(first, last)) = input1.copy();
				}
				input = af::flat(input);
			}
			else
				status = filtering(w_vec.filter, input, proj, inputScalars.Nf);
			if (status != 0)
				return -1;
		}
		if (w_vec.precondTypeMeas[0]) {
			if (DEBUG) {
				mexPrintBase("w_vec.M[timestep][subIter].dims(0) = %d\n", w_vec.M[timestep][subIter].dims(0));
				mexPrintBase("input.dims(0) = %d\n", input.dims(0));
				mexEval();
			}
			if (inputScalars.verbose >= 3)
				mexPrint("Applying diagonal normalization preconditioner (1 / (A1)), type 0");
			input /= w_vec.M[timestep][subIter];
		}
		input.eval();
		if (inputScalars.verbose >= 3)
			mexPrint("Measurement-based preconditioning applied");
	}
	return 0;
}

inline void computeIntegralImage(const scalarStruct& inputScalars, const Weighting& w_vec, const int64_t length, af::array& outputFP, af::array& meanBP) {
	if (inputScalars.BPType == 5) {
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrint("Computing integral image for backprojection");
		if (DEBUG) {
			mexPrintBase("outputFP.dims(0) = %d\n", outputFP.dims(0));
			mexPrintBase("inputScalars.nRowsD = %d\n", inputScalars.nRowsD);
			mexPrintBase("inputScalars.nColsD = %d\n", inputScalars.nColsD);
			mexPrintBase("length = %d\n", length);
			mexEval();
		}
		outputFP = af::moddims(outputFP, inputScalars.nRowsD, inputScalars.nColsD, length);
		if (inputScalars.meanBP) {
			meanBP = af::mean(af::mean(outputFP, 0), 1);
			outputFP -= af::tile(meanBP, inputScalars.nRowsD, inputScalars.nColsD, 1);
			outputFP.eval();
		}
		outputFP = af::sat(outputFP);
		outputFP = af::join(0, af::constant(0.f, 1, outputFP.dims(1), outputFP.dims(2)), outputFP);
		outputFP = af::flat(af::join(1, af::constant(0.f, outputFP.dims(0), 1, outputFP.dims(2)), outputFP));
		if (inputScalars.verbose >= 3 || DEBUG)
			mexPrint("Integral images computed");
	}
}

// Compute deconvolution (deblur/sharpening)
inline void deblur(af::array& vec, const af::array& g, const scalarStruct& inputScalars, const Weighting& w_vec, const int ii = 0) {
	if (inputScalars.verbose >= 3)
		mexPrint("Starting deconvolution");
	af::array jelppi = padding(vec, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii], inputScalars.g_dim_x + 1, inputScalars.g_dim_y + 1, inputScalars.g_dim_z + 1);
	af::array apu = padding(vec, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii], inputScalars.g_dim_x + 1, inputScalars.g_dim_y + 1, inputScalars.g_dim_z + 1);
	for (uint32_t kk = 0U; kk < inputScalars.deblur_iterations; kk++) {
		af::array apu2 = af::convolve3(jelppi, g) + inputScalars.epps;
		apu2 = apu2(af::seq(inputScalars.g_dim_x + 1, inputScalars.Nx[ii] + inputScalars.g_dim_x), af::seq(inputScalars.g_dim_y + 1, inputScalars.Ny[ii] + inputScalars.g_dim_y),
			af::seq(inputScalars.g_dim_z + 1, inputScalars.Nz[ii] + inputScalars.g_dim_z));
		apu2 = padding(apu2, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii], inputScalars.g_dim_x + 1, inputScalars.g_dim_y + 1, inputScalars.g_dim_z + 1);
		jelppi *= af::convolve3(apu / apu2, g);
		jelppi = jelppi(af::seq(inputScalars.g_dim_x + 1, inputScalars.Nx[ii] + inputScalars.g_dim_x), af::seq(inputScalars.g_dim_y + 1, inputScalars.Ny[ii] + inputScalars.g_dim_y),
			af::seq(inputScalars.g_dim_z + 1, inputScalars.Nz[ii] + inputScalars.g_dim_z));
		if (kk < inputScalars.deblur_iterations - 1)
			jelppi = padding(jelppi, inputScalars.Nx[ii], inputScalars.Ny[ii], inputScalars.Nz[ii], inputScalars.g_dim_x + 1, inputScalars.g_dim_y + 1, inputScalars.g_dim_z + 1);
	}
	vec = af::flat(jelppi);
	if (inputScalars.verbose >= 3)
		mexPrint("Deconvolution computed");
}

// The initialization steps for LSQR, CGLS, CP and FISTA algorithms
// Apply PSF blurring if applicable
inline int initializationStep(Weighting& w_vec, af::array& mData, AF_im_vectors& vec, ProjectorClass& proj, scalarStruct& inputScalars, std::vector<int64_t>& length, uint64_t m_size, const RecMethods& MethodList, uint32_t curIter, af::array& meanBP, const int64_t* pituus, const uint32_t timestep, const af::array& g = af::constant(0.f, 1, 1), const uint32_t subIter = 0, const int ii = 0, const float* atten = nullptr) {
	if (MethodList.FISTA || MethodList.FISTAL1) {
		if (curIter == 0 && subIter == 0) {
			vec.uFISTA[timestep].emplace_back(vec.im_os[timestep][ii]);
			proj.memSize += (sizeof(float) * inputScalars.im_dim[ii]);
		}
		else {
			if (inputScalars.subsetsUsed == 1 || (subIter == 0 && curIter > 0))
				vec.im_os[timestep][ii] = vec.uFISTA[timestep][ii].copy();
		}
		vec.uFISTA[timestep][ii].eval();
	}
	if (curIter == 0) {
		if (DEBUG || inputScalars.verbose >= 3)
			mexPrint("Starting initialization step");
		int status = 0;
		uint64_t yy = 0u;
		af::array mDataApu;
		if (MethodList.LSQR && subIter == 0) {
			mData = mData.as(f32);
			if (DEBUG || inputScalars.verbose >= 3)
				mexPrint("Initializing LSQR");
			vec.fLSQR[timestep].emplace_back(vec.im_os[timestep][ii].copy());
			proj.memSize += (sizeof(float) * inputScalars.im_dim[ii]);
			if (ii == 0) {
				w_vec.betaLSQR[timestep] = af::norm(mData);
				mData = mData / w_vec.betaLSQR[timestep];
			}
			if (DEBUG) {
				mexPrintBase("!!!!!!!!!!!!!!!mData = %f\n", af::sum<float>(mData));
				mexPrintBase("w_vec.betaLSQR[timestep] = %f\n", w_vec.betaLSQR[timestep]);
				mexEval();
			}
			if (inputScalars.projector_type == 6)
				backprojectionType6(mData, w_vec, vec, inputScalars, length[0], 0, proj, timestep, 0, 0, 0, 0, ii, atten);
			else {
				if (inputScalars.BPType == 5) {
					mDataApu = mData;
					computeIntegralImage(inputScalars, w_vec, length[0], mData, meanBP);
				}
				status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, MethodList, mData, 0, timestep, length, m_size, meanBP, g, proj, false, ii, pituus);
				if (status != 0) {
					return -1;
				}
				if (inputScalars.BPType == 5)
					mData = mDataApu;
			}
			if (DEBUG) {
				mexPrintBase("!!!!!!!!!!!!!!!!!!!!!!!vec.rhs_os[timestep] = %f\n", af::sum<float>(vec.rhs_os[timestep][ii]));
				mexEval();
			}
			if (ii == inputScalars.nMultiVolumes) {
				af::array temp;
				temp = vec.rhs_os[timestep][0];
				for (int ll = 1; ll <= inputScalars.nMultiVolumes; ll++)
					temp = af::join(0, temp, vec.rhs_os[timestep][ii]);
				w_vec.alphaLSQR[timestep] = af::norm(temp);
				for (int ll = 0; ll <= inputScalars.nMultiVolumes; ll++) {
					vec.im_os[timestep][ll] = vec.rhs_os[timestep][ll] / w_vec.alphaLSQR[timestep];
					vec.wLSQR[timestep].emplace_back(vec.im_os[timestep][ii].copy());
					proj.memSize += (sizeof(float) * inputScalars.im_dim[ii]);
				}
				if (DEBUG) {
					mexPrintBase("!!!!!!vec.im_os[timestep] = %f\n", af::sum<float>(vec.im_os[timestep][ii]));
					mexPrintBase("w_vec.alphaLSQR[timestep] = %f\n", w_vec.alphaLSQR[timestep]);
					mexEval();
				}
				w_vec.phiLSQR[timestep] = w_vec.betaLSQR[timestep];
				w_vec.rhoLSQR[timestep] = w_vec.alphaLSQR[timestep];
				if (DEBUG || inputScalars.verbose >= 3)
					mexPrint("LSQR initialization complete");
			}
		}
		else if (MethodList.BB && subIter == 0){
			mData = mData.as(f32);
			if (DEBUG || inputScalars.verbose >= 3)
				mexPrint("Initializing BB");
			if (inputScalars.projector_type == 6)
				backprojectionType6(mData, w_vec, vec, inputScalars, length[0], 0, proj, timestep, 0, 0, 0, 0, ii, atten);
			else {
				if (inputScalars.BPType == 5) {
					computeIntegralImage(inputScalars, w_vec, length[0], mData, meanBP);
				}
				status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, MethodList, mData, 0, timestep, length, m_size, meanBP, g, proj, false, ii, pituus);
				if (status != 0) {
					return -1;
				}
			}
			if (vec.gradBB[timestep].size() < ii + 1) {
				vec.gradBB[timestep].emplace_back(-vec.rhs_os[timestep][ii].copy());
				proj.memSize += (sizeof(float) * inputScalars.im_dim[ii]);
			}
			else
				vec.gradBB[timestep][ii] = -vec.rhs_os[timestep][ii].copy();
			if (vec.imBB[timestep].size() < ii + 1) {
				vec.imBB[timestep].emplace_back(vec.im_os[timestep][ii].copy());
				proj.memSize += (sizeof(float) * inputScalars.im_dim[ii]);
			}
			else
				vec.imBB[timestep][ii] = vec.im_os[timestep][ii].copy();
			if (w_vec.alphaBB[timestep].size() < ii + 1)
				w_vec.alphaBB[timestep].emplace_back(1e-4f);
			vec.im_os[timestep][ii] = vec.im_os[timestep][ii] - w_vec.alphaBB[timestep][ii] * vec.gradBB[timestep][ii];
			vec.gradBB[timestep][ii].eval();
			vec.im_os[timestep][ii].eval();
			if (DEBUG || inputScalars.verbose >= 3)
				mexPrint("BB initialization complete");
		}
		else if (MethodList.CGLS && subIter == 0 && !inputScalars.largeDim) {
			if (DEBUG || inputScalars.verbose >= 3)
				mexPrint("Initializing CGLS");
			mData = mData.as(f32);
			if (ii == 0)
				vec.rCGLS = mData;
			mDataApu = mData.copy();
			vec.fCGLS.emplace_back(vec.im_os[timestep][ii].copy());
			proj.memSize += (sizeof(float) * inputScalars.im_dim[ii]);
			if (inputScalars.projector_type == 6)
				backprojectionType6(mDataApu, w_vec, vec, inputScalars, length[0], 0, proj, timestep, 0, 0, 0, 0, ii, atten);
			else {
				if (inputScalars.BPType == 5) {
					computeIntegralImage(inputScalars, w_vec, length[0], mDataApu, meanBP);
				}
				status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, MethodList, mDataApu, 0, timestep, length, m_size, meanBP, g, proj, false, ii, pituus);
				if (status != 0) {
					return -1;
				}
			}
			vec.im_os[timestep][ii] = vec.rhs_os[timestep][ii].copy();
			if (ii == inputScalars.nMultiVolumes) {
				for (int ll = 0; ll <= inputScalars.nMultiVolumes; ll++)
					w_vec.gammaCGLS += af::dot<float>(vec.rhs_os[timestep][ll], vec.rhs_os[timestep][ll]);
					//w_vec.gammaCGLS += af::sum<float>(vec.rhs_os[ll] * vec.rhs_os[ll]);
				if (DEBUG || inputScalars.verbose >= 3)
					mexPrint("CGLS initialization complete");
			}
		}
		if (MethodList.SAGA && inputScalars.currentSubset == 0) {
			if (ii == 0)
				vec.stochasticHelper.resize(inputScalars.nMultiVolumes + 1);
			vec.SAGASum.emplace_back(af::constant(0.f, vec.im_os[timestep][ii].elements()));
			proj.memSize += (sizeof(float) * inputScalars.im_dim[ii]);
			for (int uu = 0; uu < inputScalars.subsetsUsed; uu++) {
				vec.stochasticHelper[ii].emplace_back(af::constant(0.f, vec.im_os[timestep][ii].elements()));
				proj.memSize += (sizeof(float) * inputScalars.im_dim[ii]);
			}
		}
		if (MethodList.CPType) {
			if (DEBUG || inputScalars.verbose >= 3)
				mexPrint("Initializing PDHG algorithm");
			if (ii == 0 && !inputScalars.largeDim && inputScalars.currentSubset == 0) {
				vec.pCP[timestep].resize(inputScalars.subsetsUsed);
				for (int uu = 0; uu < inputScalars.subsetsUsed; uu++) {
					int indD = uu + timestep * inputScalars.subsets;
					uint64_t mSize = length[indD];
					if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
						mSize = static_cast<uint64_t>(inputScalars.nRowsD) * static_cast<uint64_t>(inputScalars.nColsD) * length[indD];
					if (inputScalars.listmode && inputScalars.TOF)
						vec.pCP[timestep][uu] = af::constant(0.f, mSize);
					else
						vec.pCP[timestep][uu] = af::constant(0.f, mSize * inputScalars.nBins);
					proj.memSize += (sizeof(float) * mSize * inputScalars.nBins);
				}
			}
			else if (ii == 0 && inputScalars.largeDim)
				vec.pCP[timestep].resize(1);
			if (DEBUG) {
				mexPrintBase("subIter = %d\n", subIter);
				mexEval();
			}
			if (inputScalars.currentSubset == 0 && !inputScalars.largeDim) {
				vec.uCP[timestep].emplace_back(vec.im_os[timestep][ii].copy());
				proj.memSize += (sizeof(float) * inputScalars.im_dim[ii]);
			}
			else if (inputScalars.currentSubset == 0 && inputScalars.largeDim)
				vec.uCP[timestep].resize(1);
			if (inputScalars.verbose >= 3)
				mexPrint("PDHG initialization complete");
		}
	}
	return 0;
}

// Compute the specific weight needed by ACOSEM
inline int computeACOSEMWeight(scalarStruct& inputScalars, std::vector<int64_t>& length, float& uu, uint32_t osa_iter, uint32_t timestep, const af::array& mData, uint64_t m_size, Weighting& w_vec, AF_im_vectors& vec,
	ProjectorClass& proj, const int64_t subSum, const int64_t* pituus, const af::array& g = af::constant(0.f, 1, 1)) { // TODO vectorize ACOSEM_rhs
	if (inputScalars.verbose >= 3)
		mexPrint("Computing ACOSEM weight");
	int status = 0;
	uu = af::sum<float>(mData.as(f32));
	af::array outputFP;
	if (inputScalars.projector_type != 6) {
		if (inputScalars.listmode == 0)
			outputFP = af::constant(0.f, m_size * inputScalars.nBins);
		else
			outputFP = af::constant(0.f, m_size);
		status = forwardProjectionAFOpenCL(vec, inputScalars, w_vec, outputFP, osa_iter, timestep, length, g, m_size, proj, 0, pituus);
		if (status != 0)
			return -1;
	}
	else {
		outputFP = af::constant(0.f, inputScalars.nRowsD, inputScalars.nColsD, length[osa_iter + inputScalars.subsets * timestep]);
		forwardProjectionType6(outputFP, w_vec, vec, inputScalars, length[osa_iter + inputScalars.subsets * timestep], subSum, proj, 0, nullptr, timestep);
	}
	if (inputScalars.CT)
		w_vec.ACOSEM_rhs = af::sum<float>(af::exp(-outputFP));
	else
		w_vec.ACOSEM_rhs = af::sum<float>(outputFP);
	if (inputScalars.verbose >= 3)
		mexPrint("ACOSEM weight computed");
	return 0;
}

// The power method
inline int powerMethod(scalarStruct& inputScalars, Weighting& w_vec, std::vector<int64_t>& length, ProjectorClass& proj,
	AF_im_vectors& vec, const RecMethods& MethodList, const int64_t* pituus, const uint32_t timestep, const af::array& g = af::constant(0.f, 1, 1),
	float* F = nullptr, float* apuD = nullptr, const float* atten = nullptr) { // TODO: vectorize LCP, LCP2
	int status = 0;
	std::vector<af::array> Summ;
	af::array meanBP;
	uint64_t m_size = length[0];
	if (inputScalars.verbose > 0) {
		mexPrintBase("Starting power method for timestep %u\n", timestep);
        mexEval();
	}
	if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
		m_size = static_cast<uint64_t>(inputScalars.nRowsD) * static_cast<uint64_t>(inputScalars.nColsD) * length[0];
	std::vector<float> tauCP(inputScalars.nMultiVolumes + 1, 0.f);
	af::randomEngine r(AF_RANDOM_ENGINE_DEFAULT, 1);
	if (!inputScalars.largeDim) {
		for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
			if (ii > 0 && ii % 2 == 0 && inputScalars.im_dim[ii] == inputScalars.im_dim[ii - 1])
				vec.im_os[timestep][ii] = vec.im_os[timestep][ii - 1];
			else
				vec.im_os[timestep][ii] = af::abs(af::randn(inputScalars.im_dim[ii], f32, r));
			vec.im_os[timestep][ii] = vec.im_os[timestep][ii] / af::norm(vec.im_os[timestep][ii]);
			proj.memSize += (sizeof(float) * inputScalars.im_dim[ii]);
			vec.im_os[timestep][ii].eval();
			if (DEBUG) {
				mexPrintBase("ii = %d\n", ii);
				mexEval();
			}
		}
	}
	if (!inputScalars.largeDim) {
		for (int kk = 0; kk < w_vec.powerIterations; kk++) {
			proj.memSize += (sizeof(float) * m_size);
			af::array outputFP;
			if (inputScalars.projector_type == 6) {
				outputFP = af::constant(0.f, inputScalars.nRowsD, inputScalars.nColsD, length[0]);
				forwardProjectionType6(outputFP, w_vec, vec, inputScalars, length[0], 0, proj, 0, atten, timestep);
				outputFP.eval();
				outputFP = af::flat(outputFP);
			}
			else {
				if (inputScalars.listmode && inputScalars.TOF)
					outputFP = af::constant(0.f, m_size);
				else
					outputFP = af::constant(0.f, m_size * inputScalars.nBins);
				status = forwardProjectionAFOpenCL(vec, inputScalars, w_vec, outputFP, 0, timestep, length, g, m_size, proj, 0, pituus);
			}
			if (status != 0)
				return -1;
			if (DEBUG) {
				mexPrintBase("m_size = %d\n", m_size);
				mexEval();
			}
			status = applyMeasPreconditioning(w_vec, inputScalars, outputFP, proj, timestep);
			if (status != 0)
				return -1;
			computeIntegralImage(inputScalars, w_vec, length[0], outputFP, meanBP);
			proj.memSize += (sizeof(float) * inputScalars.im_dim[0]);
			if (inputScalars.projector_type == 6)
				backprojectionType6(outputFP, w_vec, vec, inputScalars, length[0], 0, proj, timestep, 0, 0, 0, 0, 0, atten);
			else
				status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, MethodList, outputFP, 0, timestep, length, m_size, meanBP, g, proj, false, 0, pituus);
			if (status != 0)
				return -1;
			status = applyImagePreconditioning(w_vec, inputScalars, vec.rhs_os[timestep][0], vec.im_os[timestep][0], proj, timestep, kk, 0);
			if (status != 0)
				return -1;
			tauCP[0] = (af::dot<float>(vec.im_os[timestep][0], vec.rhs_os[timestep][0]) * static_cast<float>(inputScalars.subsets)) / (af::dot<float>(vec.im_os[timestep][0], vec.im_os[timestep][0]));
			vec.im_os[timestep][0] = vec.rhs_os[timestep][0].copy();
			vec.im_os[timestep][0] /= af::norm(vec.im_os[timestep][0]);
			vec.im_os[timestep][0].eval();
			vec.rhs_os[timestep][0].eval();
			if ((inputScalars.verbose >= 2 || DEBUG)) {
				mexPrintBase("Largest eigenvalue for the main volume at iteration %d is %f\n", kk, tauCP[0]);
				mexEval();
			}
			proj.memSize -= (sizeof(float) * inputScalars.im_dim[0]);
			proj.memSize -= (sizeof(float) * m_size);
		}
		if (inputScalars.nMultiVolumes > 0) {
			for (int kk = 0; kk < w_vec.powerIterations; kk++) {
				proj.memSize += (sizeof(float) * m_size);
				af::array outputFP;
				if (inputScalars.projector_type == 6)
					outputFP = af::constant(0.f, inputScalars.nRowsD, inputScalars.nColsD, length[0]);
				else
					if (inputScalars.listmode && inputScalars.TOF)
						outputFP = af::constant(0.f, m_size);
					else
						outputFP = af::constant(0.f, m_size * inputScalars.nBins);
				if (inputScalars.projector_type == 6) {
					// Every volume projects into the same measurement, so the result is flattened
					// only after all of them are done
					forwardProjectionType6All(outputFP, w_vec, vec, inputScalars, length[0], 0, proj, atten, timestep);
				}
				else {
					for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
						status = forwardProjectionAFOpenCL(vec, inputScalars, w_vec, outputFP, 0, timestep, length, g, m_size, proj, ii, pituus);
						if (status != 0)
							return -1;
					}
				}
				if (inputScalars.projector_type == 6) {
					outputFP.eval();
					outputFP = af::flat(outputFP);
				}
				if (DEBUG) {
					mexPrintBase("m_size = %d\n", m_size);
					mexEval();
				}
				status = applyMeasPreconditioning(w_vec, inputScalars, outputFP, proj, timestep);
				if (status != 0)
					return -1;
				computeIntegralImage(inputScalars, w_vec, length[0], outputFP, meanBP);
				for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
					proj.memSize += (sizeof(float) * inputScalars.im_dim[ii]);
					if (inputScalars.projector_type == 6)
						backprojectionType6(outputFP, w_vec, vec, inputScalars, length[0], 0, proj, timestep, 0, 0, 0, 0, ii, atten);
					else
						status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, MethodList, outputFP, 0, timestep, length, m_size, meanBP, g, proj, false, ii, pituus);
					if (status != 0)
						return -1;
					if (ii == 0) {
						status = applyImagePreconditioning(w_vec, inputScalars, vec.rhs_os[timestep][ii], vec.im_os[timestep][ii], proj, timestep, kk, ii);
						if (status != 0)
							return -1;
					}
					if (ii > 0)
						tauCP[ii] = (af::dot<float>(vec.im_os[timestep][ii], vec.rhs_os[timestep][ii]) * static_cast<float>(inputScalars.subsets)) / (af::dot<float>(vec.im_os[timestep][ii], vec.im_os[timestep][ii]));
					vec.im_os[timestep][ii] = vec.rhs_os[timestep][ii].copy();
					vec.im_os[timestep][ii] /= af::norm(vec.im_os[timestep][ii]);
					vec.im_os[timestep][ii].eval();
					vec.rhs_os[timestep][ii].eval();
					if ((inputScalars.verbose >= 2 || DEBUG) && ii > 0) {
						mexPrintBase("Largest eigenvalue for volume %d at iteration %d is %f\n", ii, kk, tauCP[ii]);
						mexEval();
					}
					proj.memSize -= (sizeof(float) * inputScalars.im_dim[ii]);
				}
				proj.memSize -= (sizeof(float) * m_size);
			}
		}
	} else {
		// Random initial vector, normalized over the whole volume (not per subvolume)
		double normSq = 0.;
		for (int ii = 0; ii < inputScalars.subsets; ii++) {
			vec.im_os[timestep][0] = af::abs(af::randn(inputScalars.lDimStruct.imDim[ii], f32, r));
			normSq += static_cast<double>(af::dot<float>(vec.im_os[timestep][0], vec.im_os[timestep][0]));
			vec.im_os[timestep][0].host(&F[inputScalars.lDimStruct.cumDim[ii]]);
		}
		float scaleF = static_cast<float>(1. / std::sqrt(normSq));
		for (int kk = 0; kk < w_vec.powerIterations; kk++) {
			tauCP[0] = 0.f;
			af::array outputFP = af::constant(0.f, m_size * inputScalars.nBins);
			if (inputScalars.listmode && inputScalars.TOF)
				outputFP = af::constant(0.f, m_size);
			if (DEBUG) {
				mexPrint("Starting largeDim\n");
			}
			for (int ii = 0; ii < inputScalars.subsets; ii++) {
				largeDimFirst(inputScalars, proj, ii);
				vec.im_os[timestep][0] = af::array(inputScalars.lDimStruct.imDim[ii], &F[inputScalars.lDimStruct.cumDim[ii]], afHost) * scaleF;
				status = forwardProjectionAFOpenCL(vec, inputScalars, w_vec, outputFP, 0, timestep, length, g, m_size, proj, 0, pituus);
				if (status != 0)
					return -1;
			}
			largeDimLast(inputScalars, proj);
			status = applyMeasPreconditioning(w_vec, inputScalars, outputFP, proj, timestep);
			if (status != 0)
				return -1;
			float upper = 0.f, lower = 0.f;
			normSq = 0.;
			for (int ii = 0; ii < inputScalars.subsets; ii++) {
				largeDimFirst(inputScalars, proj, ii);
				vec.im_os[timestep][0] = af::array(inputScalars.lDimStruct.imDim[ii], &F[inputScalars.lDimStruct.cumDim[ii]], afHost) * scaleF;
				status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, MethodList, outputFP, 0, timestep, length, m_size, meanBP, g, proj, false, 0, pituus);
				if (status != 0)
					return -1;
				if (w_vec.computeD)
					w_vec.D[0][0] = af::array(inputScalars.lDimStruct.imDim[ii], &apuD[inputScalars.lDimStruct.cumDim[ii]], afHost);
				status = applyImagePreconditioning(w_vec, inputScalars, vec.rhs_os[timestep][0], vec.im_os[timestep][0], proj, timestep, kk, 0);
				if (status != 0)
					return -1;
				upper += af::dot<float>(vec.im_os[timestep][0], vec.rhs_os[timestep][0]);
				lower += af::dot<float>(vec.im_os[timestep][0], vec.im_os[timestep][0]);
				normSq += static_cast<double>(af::dot<float>(vec.rhs_os[timestep][0], vec.rhs_os[timestep][0]));
				vec.rhs_os[timestep][0].host(&F[inputScalars.lDimStruct.cumDim[ii]]);
			}
			// The stored iterate is unnormalized; it is normalized globally (scaleF) when loaded in the next iteration
			scaleF = static_cast<float>(1. / std::sqrt(normSq));
			tauCP[0] = upper * static_cast<float>(inputScalars.subsets) / lower;
			if (inputScalars.verbose >= 2 || DEBUG) {
				if (w_vec.filterIter > 0 && (w_vec.precondTypeMeas[1] || w_vec.precondTypeIm[5]))
					mexPrintBase("Largest eigenvalue at iteration %d with filtering is %f\n", kk, tauCP[0]);
				else
					mexPrintBase("Largest eigenvalue at iteration %d is %f\n", kk, tauCP[0]);
				mexEval();
			}
		}
		if (inputScalars.nMultiVolumes > 0) {
			for (int ii = 1; ii <= inputScalars.nMultiVolumes; ii++) {
				if (ii > 0 && ii % 2 == 0)
					vec.im_os[timestep][ii] = vec.im_os[timestep][ii - 1];
				else
					vec.im_os[timestep][ii] = af::abs(af::randn(inputScalars.im_dim[ii], f32, r));
				vec.im_os[timestep][ii] = vec.im_os[timestep][ii] / af::norm(vec.im_os[timestep][ii]);
				proj.memSize += (sizeof(float) * inputScalars.im_dim[ii]);
				vec.im_os[timestep][ii].eval();
				if (DEBUG) {
					mexPrintBase("ii = %d\n", ii);
					mexEval();
				}
			}
			for (int kk = 0; kk < w_vec.powerIterations; kk++) {
				proj.memSize += (sizeof(float) * m_size);
				af::array outputFP;
				if (inputScalars.projector_type == 6)
					outputFP = af::constant(0.f, inputScalars.nRowsD, inputScalars.nColsD, length[0]);
				else
					if (inputScalars.listmode && inputScalars.TOF)
						outputFP = af::constant(0.f, m_size);
					else
						outputFP = af::constant(0.f, m_size * inputScalars.nBins);
				for (int ii = 1; ii <= inputScalars.nMultiVolumes; ii++) {
					if (inputScalars.projector_type == 6) {
						// Every volume projects into the same measurement, so the result is flattened
						// only after all of them are done
						forwardProjectionType6(outputFP, w_vec, vec, inputScalars, length[0], 0, proj, ii, atten, timestep);
					}
					else {
						status = forwardProjectionAFOpenCL(vec, inputScalars, w_vec, outputFP, 0, timestep, length, g, m_size, proj, ii, pituus);
					}
					if (status != 0)
						return -1;
				}
				if (inputScalars.projector_type == 6) {
					outputFP.eval();
					outputFP = af::flat(outputFP);
				}
				if (DEBUG) {
					mexPrintBase("m_size = %d\n", m_size);
					mexEval();
				}
				status = applyMeasPreconditioning(w_vec, inputScalars, outputFP, proj, timestep);
				if (status != 0)
					return -1;
				computeIntegralImage(inputScalars, w_vec, length[0], outputFP, meanBP);
				for (int ii = 1; ii <= inputScalars.nMultiVolumes; ii++) {
					proj.memSize += (sizeof(float) * inputScalars.im_dim[ii]);
					if (inputScalars.projector_type == 6)
						backprojectionType6(outputFP, w_vec, vec, inputScalars, length[0], 0, proj, timestep, 0, 0, 0, 0, ii, atten);
					else
						status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, MethodList, outputFP, 0, timestep, length, m_size, meanBP, g, proj, false, ii, pituus);
					if (status != 0)
						return -1;
					if (ii == 0) {
						status = applyImagePreconditioning(w_vec, inputScalars, vec.rhs_os[timestep][ii], vec.im_os[timestep][ii], proj, timestep, kk, ii);
						if (status != 0)
							return -1;
					}
					if (ii > 0)
						tauCP[ii] = (af::dot<float>(vec.im_os[timestep][ii], vec.rhs_os[timestep][ii]) * static_cast<float>(inputScalars.subsets)) / (af::dot<float>(vec.im_os[timestep][ii], vec.im_os[timestep][ii]));
					vec.im_os[timestep][ii] = vec.rhs_os[timestep][ii].copy();
					vec.im_os[timestep][ii] /= af::norm(vec.im_os[timestep][ii]);
					vec.im_os[timestep][ii].eval();
					vec.rhs_os[timestep][ii].eval();
					if ((inputScalars.verbose >= 2 || DEBUG) && ii > 0) {
						mexPrintBase("Largest eigenvalue for volume %d at iteration %d is %f\n", ii, kk, tauCP[ii]);
						mexEval();
					}
					proj.memSize -= (sizeof(float) * inputScalars.im_dim[ii]);
				}
				proj.memSize -= (sizeof(float) * m_size);
			}
		}
	}
	w_vec.tauCP[timestep] = tauCP;
	uint32_t subsets = 1;
	if (!inputScalars.stochastic)
		subsets = inputScalars.subsets;
	if (w_vec.filterIter > 0 && (w_vec.precondTypeMeas[1] || w_vec.precondTypeIm[5]) && w_vec.filterIter < subsets * inputScalars.Niter) {
		bool apuM = false;
		bool apuI = false;
		if (w_vec.precondTypeMeas[1]) {
			apuM = w_vec.precondTypeMeas[1];
			w_vec.precondTypeMeas[1] = false;
		}
		else if (w_vec.precondTypeIm[5]) {
			apuI = w_vec.precondTypeIm[5];
			w_vec.precondTypeIm[5] = false;
		}
		if (inputScalars.largeDim) {
			for (int kk = 0; kk < w_vec.powerIterations; kk++) {
				tauCP[0] = 0.f;
				af::array outputFP = af::constant(0.f, m_size * inputScalars.nBins);
				if (inputScalars.listmode && inputScalars.TOF)
					outputFP = af::constant(0.f, m_size);
				if (DEBUG) {
					mexPrint("Starting largeDim\n");
				}
				for (int ii = 0; ii < inputScalars.subsets; ii++) {
					largeDimFirst(inputScalars, proj, ii);
					if (kk == 0) {
						vec.im_os[timestep][0] = af::abs(af::randn(inputScalars.lDimStruct.imDim[ii], f32, r));
						vec.im_os[timestep][0] = vec.im_os[timestep][0] / (af::norm(vec.im_os[timestep][0]) * static_cast<float>(inputScalars.subsets));
						vec.im_os[timestep][0].host(&F[inputScalars.lDimStruct.cumDim[ii]]);
					}
					else
						vec.im_os[timestep][0] = af::array(inputScalars.lDimStruct.imDim[ii], &F[inputScalars.lDimStruct.cumDim[ii]], afHost);
					status = forwardProjectionAFOpenCL(vec, inputScalars, w_vec, outputFP, 0, timestep, length, g, m_size, proj, 0, pituus);
					if (status != 0)
						return -1;
				}
				largeDimLast(inputScalars, proj);
				if (status != 0)
					return -1;
				float upper = 0.f, lower = 0.f;
				for (int ii = 0; ii < inputScalars.subsets; ii++) {
					largeDimFirst(inputScalars, proj, ii);
					vec.im_os[timestep][0] = af::array(inputScalars.lDimStruct.imDim[ii], &F[inputScalars.lDimStruct.cumDim[ii]], afHost);
					status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, MethodList, outputFP, 0, timestep, length, m_size, meanBP, g, proj, false, 0, pituus);
					if (status != 0)
						return -1;
					if (w_vec.computeD)
						w_vec.D[0][0] = af::array(inputScalars.lDimStruct.imDim[ii], &apuD[inputScalars.lDimStruct.cumDim[ii]], afHost);
					status = applyImagePreconditioning(w_vec, inputScalars, vec.rhs_os[timestep][0], vec.im_os[timestep][0], proj, timestep, kk, 0);
					if (status != 0)
						return -1;
					upper += af::dot<float>(vec.im_os[timestep][0], vec.rhs_os[timestep][0]);
					lower += af::dot<float>(vec.im_os[timestep][0], vec.im_os[timestep][0]);
					vec.im_os[timestep][0] = vec.rhs_os[timestep][0].copy();
					vec.im_os[timestep][0] /= (af::norm(vec.im_os[timestep][0]) * static_cast<float>(inputScalars.subsets));
					vec.im_os[timestep][0].host(&F[inputScalars.lDimStruct.cumDim[ii]]);
					vec.im_os[timestep][0].eval();
					vec.rhs_os[timestep][0].eval();
				}
				tauCP[0] = upper * static_cast<float>(inputScalars.subsets) / lower;
				if (inputScalars.verbose >= 2 || DEBUG) {
					mexPrintBase("Largest eigenvalue for the main volume without filtering at iteration %d is %f\n", kk, tauCP[0]);
					mexEval();
				}
			}
		}
		else {
			for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
				if (ii > 0 && ii % 2 == 0)
					vec.im_os[timestep][ii] = vec.im_os[timestep][ii - 1];
				else
					vec.im_os[timestep][ii] = af::abs(af::randn(inputScalars.im_dim[ii], f32, r));
				vec.im_os[timestep][ii] = vec.im_os[timestep][ii] / af::norm(vec.im_os[timestep][ii]);
			}
			for (int kk = 0; kk < w_vec.powerIterations; kk++) {
				af::array outputFP;
				if (inputScalars.projector_type == 6) {
					outputFP = af::constant(0.f, inputScalars.nRowsD, inputScalars.nColsD, length[0]);
					forwardProjectionType6(outputFP, w_vec, vec, inputScalars, length[0], 0, proj, 0, atten, timestep);
					outputFP.eval();
					outputFP = af::flat(outputFP);
				}
				else {
					if (inputScalars.listmode && inputScalars.TOF)
						outputFP = af::constant(0.f, m_size);
					else
						outputFP = af::constant(0.f, m_size * inputScalars.nBins);
					status = forwardProjectionAFOpenCL(vec, inputScalars, w_vec, outputFP, 0, timestep, length, g, m_size, proj, 0, pituus);
				}
				if (status != 0)
					return -1;
				if (DEBUG) {
					mexPrint("Power forward projection complete\n");
					mexEval();
				}
				status = applyMeasPreconditioning(w_vec, inputScalars, outputFP, proj, timestep);
				if (status != 0)
					return -1;
				computeIntegralImage(inputScalars, w_vec, length[0], outputFP, meanBP);
				if (inputScalars.projector_type == 6)
					backprojectionType6(outputFP, w_vec, vec, inputScalars, length[0], 0, proj, timestep, 0, 0, 0, 0, 0, atten);
				else
					status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, MethodList, outputFP, 0, timestep, length, m_size, meanBP, g, proj, false, 0, pituus);
				if (status != 0)
					return -1;
				status = applyImagePreconditioning(w_vec, inputScalars, vec.rhs_os[timestep][0], vec.im_os[timestep][0], proj, timestep, kk, 0);
				if (status != 0)
					return -1;
				tauCP[0] = (af::dot<float>(vec.im_os[timestep][0], vec.rhs_os[timestep][0]) * static_cast<float>(inputScalars.subsetsUsed)) / (af::dot<float>(vec.im_os[timestep][0], vec.im_os[timestep][0]));
				vec.im_os[timestep][0] = vec.rhs_os[timestep][0].copy();
				vec.im_os[timestep][0] /= af::norm(vec.im_os[timestep][0]);
				vec.im_os[timestep][0].eval();
				vec.rhs_os[timestep][0].eval();
				if (inputScalars.verbose >= 2 || DEBUG) {
					mexPrintBase("Largest eigenvalue for the main volume without filtering at iteration %d is %f\n", kk, tauCP[0]);
					mexEval();
				}
			}
		}
		if (inputScalars.nMultiVolumes > 0) {
			for (int kk = 0; kk < w_vec.powerIterations; kk++) {
				af::array outputFP;
				if (inputScalars.projector_type == 6)
					outputFP = af::constant(0.f, inputScalars.nRowsD, inputScalars.nColsD, length[0]);
				else
					if (inputScalars.listmode && inputScalars.TOF)
						outputFP = af::constant(0.f, m_size);
					else
						outputFP = af::constant(0.f, m_size * inputScalars.nBins);
				for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
					if (inputScalars.projector_type == 6) {
						// Every volume projects into the same measurement, so the result is flattened
						// only after all of them are done
						forwardProjectionType6(outputFP, w_vec, vec, inputScalars, length[0], 0, proj, ii, atten, timestep);
					}
					else {
						status = forwardProjectionAFOpenCL(vec, inputScalars, w_vec, outputFP, 0, timestep, length, g, m_size, proj, ii, pituus);
					}
					if (status != 0)
						return -1;
					if (DEBUG) {
						mexPrint("Power forward projection complete\n");
						mexEval();
					}
				}
				if (inputScalars.projector_type == 6) {
					outputFP.eval();
					outputFP = af::flat(outputFP);
				}
				status = applyMeasPreconditioning(w_vec, inputScalars, outputFP, proj, timestep);
				if (status != 0)
					return -1;
				computeIntegralImage(inputScalars, w_vec, length[0], outputFP, meanBP);
				for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
					if (inputScalars.projector_type == 6)
						backprojectionType6(outputFP, w_vec, vec, inputScalars, length[0], 0, proj, timestep, 0, 0, 0, 0, ii, atten);
					else
						status = backwardProjectionAFOpenCL(vec, inputScalars, w_vec, MethodList, outputFP, 0, timestep, length, m_size, meanBP, g, proj, false, ii, pituus);
					if (status != 0)
						return -1;
					if (ii == 0) {
						status = applyImagePreconditioning(w_vec, inputScalars, vec.rhs_os[timestep][ii], vec.im_os[timestep][ii], proj, timestep, kk, ii);
						if (status != 0)
							return -1;
					}
					if (ii > 0)
						tauCP[ii] = (af::dot<float>(vec.im_os[timestep][ii], vec.rhs_os[timestep][ii]) * static_cast<float>(inputScalars.subsets)) / (af::dot<float>(vec.im_os[timestep][ii], vec.im_os[timestep][ii]));
					vec.im_os[timestep][ii] = vec.rhs_os[timestep][ii].copy();
					vec.im_os[timestep][ii] /= af::norm(vec.im_os[timestep][ii]);
					vec.im_os[timestep][ii].eval();
					vec.rhs_os[timestep][ii].eval();
					if (inputScalars.verbose >= 2 || DEBUG) {
						mexPrintBase("Largest eigenvalue for volume %d without filtering at iteration %d is %f\n", ii, kk, tauCP[ii]);
						mexEval();
					}
				}
			}
		}
        w_vec.tauCP2[timestep] = tauCP;
		if (apuM) {
			w_vec.precondTypeMeas[1] = apuM;
		}
		else if (apuI) {
			w_vec.precondTypeIm[5] = apuI;
		}
	}
	for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
		proj.memSize -= (sizeof(float) * inputScalars.im_dim[ii]);
		w_vec.sigmaCP[timestep][ii] = 1.f;
		if (ii > 0)
			w_vec.sigma2CP[timestep][ii] = 1.f;
		if (inputScalars.verbose > 0) {
			if (ii == 0) {
				if (w_vec.filterIter > 0 && (w_vec.precondTypeMeas[1] || w_vec.precondTypeIm[5]))
					mexPrintBase("Largest eigenvalue for main volume with filtering is %f\n", w_vec.tauCP[timestep][ii]);
				else
					mexPrintBase("Largest eigenvalue for main volume is %f\n", w_vec.tauCP[timestep][ii]);
				if (w_vec.filterIter > 0 && (w_vec.precondTypeMeas[1] || w_vec.precondTypeIm[5]) && w_vec.filterIter < inputScalars.subsets * inputScalars.Niter)
					mexPrintBase("Largest eigenvalue for main volume without filtering is %f\n", w_vec.tauCP2[timestep][ii]);
			}
			else {
				if (w_vec.filterIter > 0 && (w_vec.precondTypeMeas[1] || w_vec.precondTypeIm[5]))
					mexPrintBase("Largest eigenvalue for volume %d with filtering is %f\n", ii, w_vec.tauCP[timestep][ii]);
				else
					mexPrintBase("Largest eigenvalue for volume %d is %f\n", ii, w_vec.tauCP[timestep][ii]);
				if (w_vec.filterIter > 0 && (w_vec.precondTypeMeas[1] || w_vec.precondTypeIm[5]) && w_vec.filterIter < inputScalars.subsets * inputScalars.Niter)
					mexPrintBase("Largest eigenvalue for volume %d without filtering is %f\n", ii, w_vec.tauCP2[timestep][ii]);
			}
			mexEval();
		}
		w_vec.LCP[timestep].emplace_back(w_vec.tauCP[timestep][ii]);
		w_vec.tauCP[timestep][ii] = 1.f / w_vec.tauCP[timestep][ii];
		if (w_vec.filterIter > 0 && (w_vec.precondTypeMeas[1] || w_vec.precondTypeIm[5])) {
			w_vec.LCP2[timestep].emplace_back(w_vec.tauCP2[timestep][ii]);
			w_vec.tauCP2[timestep][ii] = 1.f / w_vec.tauCP2[timestep][ii];
		}
	}
	return 0;
}

// Allocate vectors for proximal TV, TGV and PDHG
inline void initializeProxPriors(const RecMethods& MethodList, const scalarStruct& inputScalars, AF_im_vectors& vec, const uint32_t timestep) {
	if (MethodList.ProxTV || MethodList.ProxTGV) {
		vec.qProxTV.resize(3);
		std::fill(vec.qProxTV.begin(), vec.qProxTV.end(), af::constant(0.f, static_cast<dim_t>(inputScalars.NxPrior) * static_cast<dim_t>(inputScalars.NyPrior) * static_cast<dim_t>(inputScalars.NzPrior)));
		for (int kk = 0; kk < vec.qProxTV.size(); kk++) {
			vec.qProxTV[kk].eval();
		}
	}
	if (MethodList.ProxRDP || MethodList.ProxNLM) {
		vec.qProx.resize(1);
		std::fill(vec.qProx.begin(), vec.qProx.end(), af::constant(0.f, static_cast<dim_t>(inputScalars.NxPrior) * static_cast<dim_t>(inputScalars.NyPrior) * static_cast<dim_t>(inputScalars.NzPrior)));
		for (int kk = 0; kk < vec.qProx.size(); kk++) {
			vec.qProx[kk].eval();
		}
	}
	if (MethodList.ProxTGV) {
		if (inputScalars.TGV2D) {
			vec.vProxTGV.resize(2);
			vec.qProxTGV.resize(3);
		}
		else {
			vec.vProxTGV.resize(3);
			vec.qProxTGV.resize(6);
		}
		std::fill(vec.qProxTGV.begin(), vec.qProxTGV.end(), af::constant(0.f, static_cast<dim_t>(inputScalars.NxPrior) * static_cast<dim_t>(inputScalars.NyPrior) * static_cast<dim_t>(inputScalars.NzPrior)));
		std::fill(vec.vProxTGV.begin(), vec.vProxTGV.end(), af::constant(0.f, static_cast<dim_t>(inputScalars.NxPrior) * static_cast<dim_t>(inputScalars.NyPrior) * static_cast<dim_t>(inputScalars.NzPrior)));
		for (int kk = 0; kk < vec.qProxTGV.size(); kk++) {
			vec.qProxTGV[kk].eval();
		}
		for (int kk = 0; kk < vec.vProxTGV.size(); kk++) {
			vec.vProxTGV[kk].eval();
		}
	}
	if (MethodList.CPType && inputScalars.adaptiveType >= 1)
		vec.rhsCP[timestep].resize(inputScalars.nMultiVolumes + 1);
}

// Never-seen voxels. The OSEM-type update im * (max(rhs, epps) + epps) / (max(Summ, epps) + epps) leaves a voxel that no subset sees (the sensitivity image is
// at the epps floor in every subset) at its current value, so such voxels would keep their initial value (e.g. x0 = 1) in the output and inflate the image total
// (e.g. the axial slices beyond the detector, or the voxels of a coarse multi-resolution volume that no view reaches). After the finalization of the sensitivity
// image (every projector clamps it with Summ < epps -> epps) a voxel is seen by a subset only if Summ is above the floor. The mask (OR over the subsets) is built
// incrementally, one boolean per voxel and volume, and is applied only to the output (the final image and the saved iterations): the iterate itself is left
// untouched, so the seen voxels are bitwise identical to a run without the mask. Always on, there is no option.
// Only the algorithms that update with the per-subset sensitivity image (OSEM, ECOSEM, ROSEM, RBI, OSL-RBI, DRAMA, OSL-OSEM, ROSEM-MAP, and MLEM = OSEM with one subset) use
// the mask, the others (BSREM, PKMA, PDHG, ...) do not use Summ in this way. CT is excluded (its "sensitivity" is the backprojection of the measurement), and so is largeDim,
// where the sensitivity image is not stored per volume
inline bool seenMaskUsed(const RecMethods& MethodList, const scalarStruct& inputScalars) {
	return !inputScalars.CT && !inputScalars.largeDim && !inputScalars.fastPDHG &&
		(MethodList.OSEM || MethodList.ECOSEM || MethodList.ROSEM || MethodList.RBI || MethodList.RBIOSL || MethodList.DRAMA || MethodList.OSLOSEM || MethodList.ROSEMMAP);
}

// Adds the finalized sensitivity image Sens of the subset osa_iter (volume ii, sensitivity slot sIdx) to the never-seen voxel mask. Every subset is accumulated only
// once (the sensitivity image does not change between the iterations), so this is a no-op after the first iteration
inline void accumulateSeenMask(AF_im_vectors& vec, const scalarStruct& inputScalars, const af::array& Sens, const uint32_t sIdx, const int ii, const uint32_t osa_iter) {
	if (Sens.type() != f32 || Sens.elements() != static_cast<dim_t>(inputScalars.im_dim[ii]))
		return;
	if (vec.seenMask.size() <= sIdx) {
		vec.seenMask.resize(sIdx + 1u);
		vec.seenDone.resize(sIdx + 1u);
	}
	if (vec.seenMask[sIdx].size() <= static_cast<size_t>(ii)) {
		vec.seenMask[sIdx].resize(ii + 1);
		vec.seenDone[sIdx].resize(ii + 1);
	}
	std::vector<char>& done = vec.seenDone[sIdx][ii];
	if (done.size() < inputScalars.subsetsUsed)
		done.resize(inputScalars.subsetsUsed, 0);
	if (osa_iter >= done.size() || done[osa_iter])
		return;
	// A small relative margin so that a sensitivity at the floor is never counted as seen
	const float threshold = inputScalars.epps * (1.f + 1e-6f);
	const af::array seenNow = af::flat(Sens) > threshold;
	if (vec.seenMask[sIdx][ii].isempty())
		vec.seenMask[sIdx][ii] = seenNow;
	else
		vec.seenMask[sIdx][ii] = vec.seenMask[sIdx][ii] || seenNow;
	vec.seenMask[sIdx][ii].eval();
	done[osa_iter] = 1;
}

// Zeroes the never-seen voxels of the estimate im (volume ii of the timestep) if the mask is complete, i.e. every subset has been accumulated
inline void applySeenMask(af::array& im, const AF_im_vectors& vec, const scalarStruct& inputScalars, const uint32_t timestep, const int ii) {
	const uint32_t sIdx = inputScalars.sensPerTimestep ? timestep : 0u;
	if (sIdx >= vec.seenMask.size() || static_cast<size_t>(ii) >= vec.seenMask[sIdx].size())
		return;
	const af::array& seen = vec.seenMask[sIdx][ii];
	const std::vector<char>& done = vec.seenDone[sIdx][ii];
	if (seen.isempty() || seen.elements() != im.elements() || done.size() < inputScalars.subsetsUsed)
		return;
	for (uint32_t kk = 0; kk < inputScalars.subsetsUsed; kk++)
		if (!done[kk])
			return;
	im = af::select(af::moddims(seen, im.dims()), im, 0.f);
	im.eval();
}

// Transfer memory control back to ArrayFire
inline void transferControl(AF_im_vectors& vec, const scalarStruct& inputScalars, const af::array& g, const Weighting& w_vec, const uint32_t timestep, const uint8_t compute_norm_matrix = 2, const uint8_t no_norm = 1,
	const uint32_t osa_iter = 0, const int ii = 0) {
	// Sensitivity image slot for this timestep: 0 for a shared sensitivity image (the common case, and
	// always for Nt == 1), or the real timestep when inputScalars.sensPerTimestep requires a separate
	// slot per timestep (see reconstructionAF.h). When the caller passed no_norm == 1 because this
	// timestep's backprojection did not (re)compute the sensitivity image (shared-image case, timestep
	// other than the one that owns it for this subset), the no_norm == 0u guards below skip finalizing
	// it again -- only the unlock() runs, releasing the lock transferSensitivityImage took to pass the
	// (untouched) image to the backprojection kernel.
	const uint32_t sIdx = inputScalars.sensPerTimestep ? timestep : 0u;
	if (compute_norm_matrix == 1u) {
		vec.Summ[sIdx][ii][0].unlock();
		if (no_norm == 0u) {
#if defined(OPENCL) || defined(METAL)
			if (inputScalars.atomic_64bit)
				vec.Summ[sIdx][ii][0] = vec.Summ[sIdx][ii][0].as(f32) / TH;
			else if (inputScalars.atomic_32bit)
				vec.Summ[sIdx][ii][0] = vec.Summ[sIdx][ii][0].as(f32) / TH32;
#endif
			if (inputScalars.use_psf) {
				vec.Summ[sIdx][ii][0] = computeConvolution(vec.Summ[sIdx][ii][0], g, inputScalars, w_vec, 1, ii);
			}
			// Prevent division by zero
			vec.Summ[sIdx][ii][0](vec.Summ[sIdx][ii][0] < inputScalars.epps) = inputScalars.epps;
			vec.Summ[sIdx][ii][0].eval();
			if (DEBUG) {
				mexPrint("Sens image steps 1 done\n");
			}
		}
	}
	else if (compute_norm_matrix == 2) {
		vec.Summ[sIdx][ii][osa_iter].unlock();
		if (no_norm == 0u) {
#if defined(OPENCL) || defined(METAL)
			if (inputScalars.atomic_64bit) {
				vec.Summ[sIdx][ii][osa_iter] = vec.Summ[sIdx][ii][osa_iter].as(f32) / TH;
			}
			else if (inputScalars.atomic_32bit) {
				vec.Summ[sIdx][ii][osa_iter] = vec.Summ[sIdx][ii][osa_iter].as(f32) / TH32;
			}
#endif
			if (inputScalars.use_psf) {
				vec.Summ[sIdx][ii][osa_iter] = computeConvolution(vec.Summ[sIdx][ii][osa_iter], g, inputScalars, w_vec, 1, ii);
			}
			vec.Summ[sIdx][ii][osa_iter](vec.Summ[sIdx][ii][osa_iter] < inputScalars.epps) = inputScalars.epps;
			vec.Summ[sIdx][ii][osa_iter].eval();
			if (DEBUG) {
				mexPrint("Sens image steps 2 done\n");
			}
			if (DEBUG) {
				mexPrintBase("inputScalars.epps = %f\n", inputScalars.epps);
				mexPrintBase("min(Summ) = %f\n", af::min<float>(vec.Summ[sIdx][ii][osa_iter]));
				mexEval();
			}
		}
	}
	if (DEBUG && inputScalars.atomic_64bit) {
		mexPrintBase("min(rhs_os) = %d\n", af::min<int64_t>(vec.rhs_os[timestep][ii]));
		mexPrintBase("inputScalars.atomic_64bit = %d\n", inputScalars.atomic_64bit);
		mexEval();
	}
}
