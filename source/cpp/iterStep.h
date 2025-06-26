#pragma once
#include "functions.hpp"
#include "algorithms.h"
#include "priors.h"
#include <cstring>

template <typename T>
int computeOSEstimatesIter(AF_im_vectors& vec, Weighting& w_vec, const RecMethods& MethodList, const scalarStruct& inputScalars, const uint32_t iter,
	ProjectorClass& proj, const af::array& g, T* output, uint32_t& ee, size_t& tt, const float* x0) {

	// Compute BSREM and ROSEMMAP updates if applicable
	// Otherwise simply save the current iterate if applicable
	uint32_t it = 0U;
	//if (inputScalars.FISTAAcceleration) {
	//	//if (iter > 0 && iter >= w_vec.filterIter / inputScalars.subsets) {
	//	if (iter > 0) {
	//		const float t = w_vec.tFISTA;
	//		w_vec.tFISTA = (1.f + std::sqrt(1.f + 4.f * w_vec.tFISTA * w_vec.tFISTA)) / 2.f;
	//		//}
	//		for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
	//			const af::array apu = vec.im_os[ii].copy();
	//			vec.im_os[ii] = vec.im_os[ii] + (t - 1.f) / w_vec.tFISTA * (vec.im_os[ii] - vec.FISTAApu[ii]);
	//			af::eval(vec.im_os[ii]);
	//			vec.FISTAApu[ii] = apu.copy();
	//		}
	//	}
	//	else {
	//		w_vec.tFISTA = (1.f + std::sqrt(1.f + 4.f * w_vec.tFISTA * w_vec.tFISTA)) / 2.f;
	//			for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
	//				if (iter == 0)
	//					vec.FISTAApu.emplace_back(vec.im_os[ii].copy());
	//				else
	//					vec.FISTAApu[ii] = vec.im_os[ii].copy();
	//			}
	//	}
	//}
	if (MethodList.BSREM || MethodList.ROSEMMAP) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing regularization for BSREM/ROSEMMAP");
		int status = 0;
		const af::array im = vec.im_os[0];
		status = applyPrior(vec, w_vec, MethodList, inputScalars, proj, w_vec.beta, iter, 0, true);
		if (status != 0)
			return -1;
		// MAP/Prior-algorithms
		// Special case for BSREM and ROSEM-MAP
		vec.im_os[0] = MAP(im, w_vec.lambda[iter], vec.im_os[0], inputScalars.epps);
		if (inputScalars.verbose >= 3)
			mexPrint("Regularization for BSREM/ROSEMMAP computed");
	}
	if (inputScalars.saveIter || (inputScalars.saveIterationsMiddle > 0 && (iter == inputScalars.Niter - 1 || inputScalars.saveNIter[ee] == iter))) {
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
		float* jelppi = getSingles(output, "solu");
#else
		float* jelppi = output;
#endif
		if (inputScalars.saveIter && iter == 0) {
			std::memcpy(&jelppi[tt], &x0[0], inputScalars.im_dim[0] * sizeof(float));
			tt += inputScalars.im_dim[0];
		}
		if (inputScalars.use_psf && inputScalars.deconvolution) {
			af::array apu = vec.im_os[0].copy();
			deblur(apu, g, inputScalars, w_vec);
			if (inputScalars.useHalf)
				apu = apu.as(f32);
			apu.host(&jelppi[tt]);
		}
		else {
			if (inputScalars.useHalf) {
				af::array tempIm = vec.im_os[0].copy().as(f32);
				tempIm.host(&jelppi[tt]);
			}
			else
				vec.im_os[0].host(&jelppi[tt]);
		}
		ee++;
		tt += inputScalars.im_dim[0];
	}
	return 0;
}