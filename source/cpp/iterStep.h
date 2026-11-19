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
		// The spatial prior gradient is recomputed only every regEveryIter-th iteration and the last computed gradient
		// (vec.dU) is used in between. This function runs once per iteration (iter from 0 to inputScalars.Niter - 1);
		// the first and last are always computed.
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
		MAP(vec.im_os[timestep][0], w_vec.lambda[timestep][iter], totalPriorGradient(vec, MethodList, timestep), inputScalars.epps);
		if (inputScalars.verbose >= 3)
			mexPrint("Regularization for BSREM/ROSEMMAP computed");
	}
	// Enforce positivity if applicable. This is done before saving so that the saved and the returned image are both clamped
	if (inputScalars.enforcePositivity && inputScalars.subsetsUsed == 1 && !MethodList.CPType) {
		for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++)
			vec.im_os[timestep][ii] = af::select(vec.im_os[timestep][ii] < inputScalars.epps, inputScalars.epps, vec.im_os[timestep][ii]);
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
		// im_os holds the Krylov vector for LSQR/CGLS except after the last iteration, so the solution iterate is saved instead
		const af::array* saveSrc = &vec.im_os[timestep][0];
		if (MethodList.LSQR && vec.fLSQR.size() > timestep && !vec.fLSQR[timestep].empty())
			saveSrc = &vec.fLSQR[timestep][0];
		else if (MethodList.CGLS && !vec.fCGLS.empty())
			saveSrc = &vec.fCGLS[0];
		if (inputScalars.use_psf && inputScalars.deconvolution) {
			af::array apu = saveSrc->copy();
			deblur(apu, g, inputScalars, w_vec);
			apu.host(&jelppi[offset]);
		}
		else
			saveSrc->host(&jelppi[offset]);
	}
	return 0;
}