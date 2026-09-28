#pragma once
#include "functions.hpp"
#include "algorithms.h"
#include "priors.h"
#include <cstring>

template <typename T>
int computeOSEstimatesIter(AF_im_vectors& vec, Weighting& w_vec, const RecMethods& MethodList, const scalarStruct& inputScalars, const uint32_t iter,
	ProjectorClass& proj, const af::array& g, T* output, const bool doSave, const uint32_t slot, const float* x0, const int timestep) {

	// Compute BSREM and ROSEMMAP updates if applicable
	// Otherwise simply save the current iterate if applicable
	uint32_t it = 0U;
	if (MethodList.BSREM || MethodList.ROSEMMAP) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing regularization for BSREM/ROSEMMAP");
		int status = 0;
		// Compute the spatial prior only every regEveryIter-th iteration. This function runs once per
		// iteration (iter from 0 to inputScalars.Niter - 1); the first and last are always computed.
		const int64_t regCounter = static_cast<int64_t>(iter);
		const int64_t regCounterMax = static_cast<int64_t>(inputScalars.Niter) - 1;
		const bool computeReg = inputScalars.regEveryIter <= 1 || regCounter == 0 || regCounter == regCounterMax
			|| (regCounter % static_cast<int64_t>(inputScalars.regEveryIter)) == 0;
		if (computeReg) {
			status = applySpatialPrior(vec, w_vec, MethodList, inputScalars, proj, w_vec.beta, timestep, iter, 0, true);
			if (status != 0)
				return -1;
		}
		// MAP/Prior-algorithms
		// Special case for BSREM and ROSEM-MAP
		MAP(vec.im_os[timestep][0], w_vec.lambda[timestep][iter], vec.dU[timestep], inputScalars.epps);
		if (inputScalars.verbose >= 3)
			mexPrint("Regularization for BSREM/ROSEMMAP computed");
	}
	if (doSave) {
		if (inputScalars.verbose >= 3)
			mexPrintVar("Saving intermediate result at iteration ", iter);
		if (DEBUG) {
			mexPrintBase("iter = %d\n", iter);
			mexPrintBase("slot = %d\n", slot);
			mexEval();
		}
#ifdef MATLAB
		float* jelppi = getSingles(output, "solu");
#else
		float* jelppi = output;
#endif
		// Output memory layout is [Nx,Ny,Nz,Nt,saves], i.e. voxel fastest, then timestep, then
		// save slot: the offset for a given (slot, timestep) pair is (slot * Nt + timestep) * im_dim[0].
		// This is computed explicitly here (rather than via a running write offset) so that it is
		// correct regardless of the order in which timesteps are processed, and identical to the
		// pre-existing single-timestep (Nt == 1) layout when Nt == 1.
		const size_t timestepOffset = static_cast<size_t>(timestep) * static_cast<size_t>(inputScalars.im_dim[0]);
		const size_t slotStride = static_cast<size_t>(inputScalars.Nt) * static_cast<size_t>(inputScalars.im_dim[0]);
		if (inputScalars.saveIter && iter == 0) {
			// Slot 0 always holds the initial value x0. x0 holds a single initial image shared by
			// all timesteps (see the Nt-loop initialization above that (re)reads x0 from offset 0
			// for every timestep), so the same source data is copied into each timestep's slot-0 region.
			std::memcpy(&jelppi[timestepOffset], &x0[0], inputScalars.im_dim[0] * sizeof(float));
		}
		const size_t offset = static_cast<size_t>(slot) * slotStride + timestepOffset;
		if (inputScalars.use_psf && inputScalars.deconvolution) {
			af::array apu = vec.im_os[timestep][0].copy();
			deblur(apu, g, inputScalars, w_vec);
			apu.host(&jelppi[offset]);
		}
		else
			vec.im_os[timestep][0].host(&jelppi[offset]);
	}
	return 0;
}