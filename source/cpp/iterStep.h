#pragma once
#include "functions.hpp"
#include "algorithms.h"
#include "priors.h"
#include <cstring>

template <typename T>
int computeOSEstimatesIter(AF_im_vectors& vec, Weighting& w_vec, const RecMethods& MethodList, const scalarStruct& inputScalars, const uint32_t iter,
	ProjectorClass& proj, const af::array& g, T* output, uint32_t& ee, size_t& tt, const float* x0, const int timestep) {

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
	const bool storeMultiResolution = inputScalars.storeMultiResolution && inputScalars.nMultiVolumes > 0;
    bool saveCurrentMultiResolution = inputScalars.saveIter;
    size_t multiResolutionSaveIndex = static_cast<size_t>(iter) + 1ULL;
    if (storeMultiResolution && !inputScalars.saveIter && inputScalars.saveIterationsMiddle > 0) {
        saveCurrentMultiResolution = false;
        for (size_t ii = 0; ii < inputScalars.saveIterationsMiddle; ii++) {
            if (inputScalars.saveNIter[ii] == iter) {
                multiResolutionSaveIndex = ii;
                saveCurrentMultiResolution = true;
                break;
            }
        }
        if (!saveCurrentMultiResolution && iter == inputScalars.Niter - 1) {
            multiResolutionSaveIndex = inputScalars.saveIterationsMiddle;
            saveCurrentMultiResolution = true;
        }
    }
    const bool saveCurrentPrimary = !storeMultiResolution && (inputScalars.saveIter ||
        (inputScalars.saveIterationsMiddle > 0 && (iter == inputScalars.Niter - 1 ||
            (ee < inputScalars.saveIterationsMiddle && inputScalars.saveNIter[ee] == iter))));
    if ((storeMultiResolution && saveCurrentMultiResolution) || saveCurrentPrimary) {
		if (inputScalars.verbose >= 3)
			mexPrintVar("Saving intermediate result at iteration ", iter);
		if (DEBUG) {
			mexPrintBase("iter = %d\n", iter);
			mexPrintBase("ee = %d\n", ee);
			if (inputScalars.saveIterationsMiddle > 0)
				mexPrintBase("inputScalars.saveNIter[ee] = %d\n", inputScalars.saveNIter[ee]);
			mexEval();
		}
#ifdef MATLAB
        if (storeMultiResolution) {
            const size_t savedFrameCount = inputScalars.saveIter ? static_cast<size_t>(inputScalars.Niter) + 1ULL : inputScalars.saveIterationsMiddle + 1ULL;
            size_t x0Offset = 0ULL;
            for (uint32_t ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
                float* outputVolume = getSingles(output, "", static_cast<mwIndex>(ii));
                const size_t timestepOffset = static_cast<size_t>(timestep) * savedFrameCount * inputScalars.im_dim[ii];
                if (inputScalars.saveIter && iter == 0)
                    std::memcpy(&outputVolume[timestepOffset], &x0[x0Offset], inputScalars.im_dim[ii] * sizeof(float));
                const size_t outputOffset = timestepOffset + static_cast<size_t>(multiResolutionSaveIndex) * inputScalars.im_dim[ii];
                if (ii == 0 && inputScalars.use_psf && inputScalars.deconvolution) {
                    af::array apu = vec.im_os[timestep][ii].copy();
                    deblur(apu, g, inputScalars, w_vec);
                    apu.host(&outputVolume[outputOffset]);
                } else {
                    vec.im_os[timestep][ii].host(&outputVolume[outputOffset]);
                }
                x0Offset += inputScalars.im_dim[ii];
            }
            if (timestep == inputScalars.Nt - 1)
                ee++;
            return 0;
        }
		float* jelppi = getSingles(output, "solu");
#else
		float* jelppi = output;
#endif
		if (inputScalars.saveIter && iter == 0) {
			std::memcpy(&jelppi[tt], &x0[0], inputScalars.im_dim[0] * sizeof(float));
			tt += inputScalars.im_dim[0];
		}
		if (inputScalars.use_psf && inputScalars.deconvolution) {
			af::array apu = vec.im_os[timestep][0].copy();
			deblur(apu, g, inputScalars, w_vec);
			apu.host(&jelppi[tt]);
		}
		else
			vec.im_os[timestep][0].host(&jelppi[tt]);
		ee++;
		tt += inputScalars.im_dim[0];
	}
	return 0;
}