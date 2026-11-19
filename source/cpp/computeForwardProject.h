#pragma once
#include "functions.hpp"

// Computes all computations using the forward projection and outputting a measurement-domain vector
inline int computeForwardStep(const RecMethods& MethodList, af::array& y, af::array& input, const int64_t length, const scalarStruct& inputScalars, Weighting& w_vec, const af::array& randomsData, AF_im_vectors& vec, ProjectorClass& proj, const uint32_t timestep, const uint32_t iter = 0, const uint32_t subIter = 0, const int ii = 0, float* residual = nullptr, const af::array& aNorm = af::array()) {
	if (DEBUG || inputScalars.verbose >= 3) {
		proj.tStartLocal = std::chrono::steady_clock::now();
	}
	//af::deviceGC();
	int status = 0;
	af::array indeksit;
	bool indS = false;
	// CT host-side normalization (detector gain): exp(-FP) becomes n * exp(-FP). aNorm is empty when not used
	const bool hostGain = inputScalars.CT && !aNorm.isempty();
	uint32_t kk = inputScalars.currentSubset + iter * inputScalars.subsets;
	if (DEBUG) {
		mexPrintBase("length = %d\n", length);
		mexPrintBase("y.dims(0) = %d\n", y.dims(0));
		mexPrintBase("input.dims(0) = %d\n", input.dims(0));
		mexEval();
	}
	// LSQR, CGLS and BB are linear least-squares methods and need a linear operator: the epps only guards the
	// divisions of the Poisson methods
	if (!inputScalars.CT && !(MethodList.LSQR || MethodList.CGLS || MethodList.BB))
		input += inputScalars.epps;
	if (w_vec.precondTypeMeas[1] && w_vec.filterIter > 0 && kk == w_vec.filterIter) {
		if (inputScalars.verbose >= 3)
			mexPrint("Filter iterations complete. Switching tau/sigma-values");
		if (MethodList.CPType || MethodList.FISTA || MethodList.FISTAL1 || MethodList.ProxTGV || MethodList.ProxTV) {
			// The switch is done for every timestep, since this function is called separately for each of them
			for (uint32_t t = 0; t < inputScalars.Nt; t++) {
				for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
					if (w_vec.sigmaCP[t][ii] == 1.f)
						w_vec.tauCP[t][ii] = w_vec.tauCP2[t][ii];
					else if (w_vec.sigmaCP[t][ii] == w_vec.tauCP[t][ii]) {
						w_vec.tauCP[t][ii] = w_vec.tauCP2[t][ii];
						w_vec.sigmaCP[t][ii] = w_vec.tauCP2[t][ii];
					}
					else
						w_vec.sigmaCP[t][ii] = w_vec.tauCP2[t][ii];
					if (inputScalars.adaptiveType == 1)
						w_vec.alphaCP[t][ii] = 1.f;
					else if (inputScalars.adaptiveType == 2)
						w_vec.alphaCP[t][ii] = .95f;
				}
			}
			w_vec.LCP = w_vec.LCP2;
		}
		if (MethodList.MRAMLA || MethodList.MBSREM || MethodList.SPS || MethodList.RAMLA || MethodList.BSREM || MethodList.ROSEM || MethodList.ROSEMMAP || MethodList.PKMA || MethodList.SAGA) {
			for (uint32_t t = 0; t < inputScalars.Nt; t++)
				w_vec.lambda[t] = w_vec.lambdaFiltered[t];
		}
		w_vec.precondTypeMeas[1] = false;
	}
	// Ahn & Fessler's modified log-likelihood applies to all LORs with zero additive term. The indices must be
	// computed before the randoms/scatter are added to the forward projection
	if ((MethodList.MBSREM || MethodList.MRAMLA || MethodList.SPS) && !inputScalars.CT) {
		if (!inputScalars.randoms_correction) {
			indeksit = y > 0 && input <= w_vec.epsilon_mramla[timestep];
			indS = af::anyTrue<bool>(indeksit);
		}
		else if (!af::allTrue<bool>(randomsData > 0.f)) {
			if (inputScalars.TOF && inputScalars.listmode == 0) {
				// A TOF-resolved additive term (randomsTOF) has a value for each TOF bin, otherwise it is repeated for the bins
				// (list-mode data has one value per event, in the TOF bin of the event, for both the data and the additive term)
				const af::array indR = inputScalars.randomsTOF ? (randomsData == 0.f) : af::tile(randomsData == 0.f, inputScalars.nBins);
				indeksit = y > 0 && indR && input <= w_vec.epsilon_mramla[timestep];
			}
			else
				indeksit = y > 0 && randomsData == 0.f && input <= w_vec.epsilon_mramla[timestep];
			indS = af::anyTrue<bool>(indeksit);
		}
	}
	if (inputScalars.randoms_correction) {
		if (!inputScalars.CT) {
			if (DEBUG) {
				mexPrintBase("input.dims(0) = %d\n", input.dims(0));
				mexPrintBase("randomsData.elements() = %d\n", randomsData.elements());
				mexPrintBase("randomsData.dims(0) = %d\n", randomsData.dims(0));
				mexPrintBase("randomsData.dims(1) = %d\n", randomsData.dims(1));
				mexEval();
			}
			if (inputScalars.verbose >= 3)
				mexPrint("Adding randoms/scatter data to forward projection");
			if (inputScalars.TOF && !inputScalars.randomsTOF && inputScalars.listmode == 0)
				input += af::tile(randomsData, inputScalars.nBins);
			else
				input += randomsData;
		}
	}
	if (MethodList.CPType && inputScalars.subsetsUsed > 1) {
		// Shallow copy: pCP[timestep][subIter] is only ever rebound later (never written in place or through a raw device pointer),
		// so p0CP keeps the old values
		vec.p0CP[timestep] = vec.pCP[timestep][subIter];
	}
	if (MethodList.ACOSEM || MethodList.OSLCOSEM > 0 || MethodList.OSEM || MethodList.COSEM || MethodList.ECOSEM ||
		MethodList.ROSEM || MethodList.OSLOSEM || MethodList.ROSEMMAP) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing (OSL-)(A)(E)(C)(R)OSEM");
		if (inputScalars.CT) {
			input = af::select(input < 0.f, 0.f, input);
			if (inputScalars.verbose >= 3)
				mexPrint("CT mode");
			if (inputScalars.verbose >= 3 && inputScalars.randoms_correction)
				mexPrint("Adding scatter data to forward projection");
			if (MethodList.OSEM || MethodList.ROSEM || MethodList.OSLOSEM || MethodList.ROSEMMAP) {
				if (hostGain)
					input = aNorm * af::exp(-input);
				else
					input = af::exp(-input);
			}
			else
				input = af::exp(-input) / y.as(f32);
		}
		else {
			if (inputScalars.verbose >= 3)
				mexPrint("PET/SPECT mode");
			input = y.as(f32) / (input);
		}
	}
	else if (MethodList.RAMLA || MethodList.BSREM || MethodList.RBI || MethodList.RBIOSL || MethodList.DRAMA) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing RAMLA/BSREM/RBI/DRAMA");
		if (inputScalars.CT) {
			input = af::select(input < 0.f, 0.f, input);
			if (inputScalars.verbose >= 3)
				mexPrint("CT mode");
			if (inputScalars.verbose >= 3 && inputScalars.randoms_correction)
				mexPrint("Adding scatter data to forward projection");
			if (inputScalars.randoms_correction) {
				const af::array expInput = (hostGain ? af::array(aNorm * af::exp(-input)) : af::array(af::exp(-input)));
				input = expInput * inputScalars.flat - (y.as(f32) * expInput) / (expInput + randomsData);
			}
			else
				input = (hostGain ? af::array(aNorm * af::exp(-input)) : af::array(af::exp(-input))) * inputScalars.flat - y.as(f32);
		}
		else {
			if (inputScalars.verbose >= 3)
				mexPrint("PET/SPECT mode");
			if (inputScalars.listmode > 0)
				input = y.as(f32) / (input);
			else
				input = y.as(f32) / (input) - 1.f;
		}
	}
	else if (MethodList.PKMA) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing PKMA");
		if (inputScalars.CT) {
			input = af::select(input < 0.f, 0.f, input);
			if (inputScalars.verbose >= 3)
				mexPrint("CT mode");
			if (inputScalars.verbose >= 3 && inputScalars.randoms_correction)
				mexPrint("Adding scatter data to forward projection");
			if (inputScalars.randoms_correction) {
				const af::array expInput = (hostGain ? af::array(aNorm * af::exp(-input)) : af::array(af::exp(-input)));
				input = (y.as(f32) * expInput) / (expInput + randomsData) - expInput * inputScalars.flat;
			}
			else
				input = y.as(f32) - (hostGain ? af::array(aNorm * af::exp(-input)) : af::array(af::exp(-input))) * inputScalars.flat;
		}
		else {
			if (inputScalars.verbose >= 3)
				mexPrint("PET/SPECT mode");
			if (inputScalars.listmode > 0)
				input = y.as(f32) / (input);
			else
				input = 1.f - y.as(f32) / (input);
		}
		status = applyMeasPreconditioning(w_vec, inputScalars, input, proj, timestep, subIter);
		if (status != 0)
			return -1;
	}
	else if (MethodList.MBSREM || MethodList.MRAMLA || MethodList.SPS) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing MBSREM/MRAMLA/SPS");
		if (inputScalars.CT) {
			input = af::select(input < 0.f, 0.f, input);
			if (inputScalars.verbose >= 3)
				mexPrint("CT mode");
			if (inputScalars.verbose >= 3 && inputScalars.randoms_correction)
				mexPrint("Adding scatter data to forward projection");
			if (inputScalars.randoms_correction) {
				const af::array expInput = (hostGain ? af::array(aNorm * af::exp(-input)) : af::array(af::exp(-input)));
				input = expInput * inputScalars.flat - (y.as(f32) * expInput) / (expInput + randomsData);
			}
			else
				input = (hostGain ? af::array(aNorm * af::exp(-input)) : af::array(af::exp(-input))) * inputScalars.flat - y.as(f32);
		}
		else {
			if (inputScalars.verbose >= 3)
				mexPrint("PET/SPECT mode");
			if (indS) {
				input = af::select(indeksit,
					y.as(f32) / w_vec.epsilon_mramla[timestep] - 1.f - (y.as(f32) / (w_vec.epsilon_mramla[timestep] * w_vec.epsilon_mramla[timestep])) * (input - w_vec.epsilon_mramla[timestep]),
					y.as(f32) / (input + inputScalars.epps) - 1.f);
				// The formulas above include the -1 term. For list-mode data it is replaced by the sensitivity image
				// (subtracted later), so remove it here
				if (inputScalars.listmode > 0)
					input += 1.f;
			}
			else
				if (inputScalars.listmode > 0)
					input = y.as(f32) / (input);
				else
					input = y.as(f32) / (input) - 1.f;
		}
		status = applyMeasPreconditioning(w_vec, inputScalars, input, proj, timestep, subIter);
		if (status != 0)
			return -1;
	}
	else if (MethodList.LSQR) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing LSQR");
		input -= w_vec.alphaLSQR[timestep] * y.as(f32);
		w_vec.betaLSQR[timestep] = af::norm(input);
		input = input / w_vec.betaLSQR[timestep];
		y = input.copy();
	}
	else if (MethodList.CGLS) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing CGLS");
		const float normi = af::norm(input);
		w_vec.alphaCGLS = w_vec.gammaCGLS / (normi * normi);
		input = vec.rCGLS - w_vec.alphaCGLS * input;
		vec.rCGLS = input.copy();
	}
	else if (MethodList.BB) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing BB");
		input -= y.as(f32);
	}
	else if (MethodList.SART || MethodList.POCS) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing SART or ASD-POCS");
		if (MethodList.POCS)
			vec.f0POCS[timestep] = vec.im_os[timestep];
		input = y.as(f32) - input;
		if (inputScalars.storeResidual) {
			//residual[kk] = af::sum<float>(af::matmulTN(input, input)) * .5;
			residual[kk] = af::norm(input);
			residual[kk] = residual[kk] * residual[kk] * .5f;
		}
		input /= w_vec.M[timestep][subIter];
	}
	else if (MethodList.PDHG || MethodList.PDDY) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing PDHG");
		if (DEBUG) {
			mexPrintBase("rdim0 = %u\n", y.dims(0));
			mexPrintBase("rdim1 = %u\n", y.dims(1));
			mexEval();
		}
		af::array res = input - y.as(f32);
		if (DEBUG) {
			mexPrintBase("res = %f\n", af::sum<float>(res));
			mexPrintBase("max(res) = %f\n", af::max<float>(res));
			mexEval();
		}
		if (inputScalars.storeResidual) {
			//af::array ressa = res.copy();
			//status = applyMeasPreconditioning(w_vec, inputScalars, ressa, proj, timestep, subIter);
			//residual[kk] = af::sum<float>(af::matmulTN(res, ressa)) * .5;
			residual[kk] = af::norm(res);
			residual[kk] = residual[kk] * residual[kk] * .5f;
		}
		status = applyMeasPreconditioning(w_vec, inputScalars, res, proj, timestep, subIter);
		if (status != 0)
			return -1;
		if (DEBUG) {
			mexPrintBase("res = %f\n", af::sum<float>(res));
			mexPrintBase("max(res) = %f\n", af::max<float>(res));
			mexEval();
		}
		if (w_vec.precondTypeMeas[1] && subIter + iter * inputScalars.subsets < w_vec.filterIter) {
			if (FINVERSE) {
				if (inputScalars.verbose >= 3)
					mexPrint("Computing inverse with circulant matrix");
				input = (vec.pCP[timestep][subIter] + w_vec.sigmaCP[timestep][ii] * res);
				if (inputScalars.subsetType == 5 || inputScalars.subsetType == 4) {
					if (inputScalars.subsetType == 4)
						input = af::moddims(input, inputScalars.nRowsD, input.elements() / inputScalars.nRowsD);
					else
						input = af::moddims(input, inputScalars.nColsD, input.elements() / inputScalars.nColsD);
				}
				else
					input = af::moddims(input, inputScalars.nRowsD, inputScalars.nColsD, input.elements() / (inputScalars.nRowsD * inputScalars.nColsD));
				if (inputScalars.adaptiveType >= 1 && ii == 0) {
					w_vec.Ffilter = af::ifft(w_vec.filter) * w_vec.sigmaCP[timestep][ii];
					w_vec.Ffilter(0) = w_vec.Ffilter(0) + 1.f;
					w_vec.Ffilter = af::real(af::fft(w_vec.Ffilter));
				}
				status = filteringInverse(w_vec.Ffilter, input, proj, inputScalars.Nf);
				if (status != 0)
					return -1;
			}
			else {
				af::array apu = af::tile(w_vec.filter, res.elements() / w_vec.filter.elements());
				if (DEBUG) {
					mexPrintBase("dim(0) = %d\n", apu.dims(0));
					mexPrintBase("dim(1) = %d\n", apu.dims(1));
					mexPrintBase("vec.pCP[timestep][subIter].dim(0) = %d\n", vec.pCP[timestep][subIter].dims(0));
					mexPrintBase("res.dim(0) = %d\n", res.dims(0));
					mexPrintBase("res.elements() = %d\n", res.elements());
					mexPrintBase("apu.elements() = %d\n", apu.elements());
					mexPrintBase("res.elements() / w_vec.filter.elements() = %d\n", res.elements() / apu.elements());
					mexEval();
				}
				input = (vec.pCP[timestep][subIter] + w_vec.sigmaCP[timestep][ii] * res);
			}
		}
		else {
			// Diagonal preconditioning (Pock and Chambolle 2011) uses Sigma = sigma / M. res has already been divided by M,
			// so the numerator contains Sigma * (Ax - y). The prox of the conjugate of ||z - y||^2 / 2 is
			// (v - Sigma y) / (1 + Sigma) elementwise, i.e. the denominator must also use Sigma. Dividing by (1 + sigma)
			// instead solved the M^-1-weighted least squares problem
			if (MethodList.ProxTGV) {
				if (inputScalars.verbose >= 3)
					mexPrint("Computing Proximal TGV");
				if (w_vec.precondTypeMeas[0])
					input = (vec.pCP[timestep][subIter] + w_vec.sigmaCP[timestep][ii] * res) / (1.f + w_vec.sigmaCP[timestep][ii] * w_vec.betaReg / w_vec.M[timestep][subIter]);
				else
					input = (vec.pCP[timestep][subIter] + w_vec.sigmaCP[timestep][ii] * res) / (1.f + w_vec.sigmaCP[timestep][ii] * w_vec.betaReg);
			}
			else {
				if (inputScalars.verbose >= 3)
					mexPrint("Computing PDHG");
				if (w_vec.precondTypeMeas[0])
					input = (vec.pCP[timestep][subIter] + w_vec.sigmaCP[timestep][ii] * res) / (1.f + w_vec.sigmaCP[timestep][ii] / w_vec.M[timestep][subIter]);
				else
					input = (vec.pCP[timestep][subIter] + w_vec.sigmaCP[timestep][ii] * res) / (1.f + w_vec.sigmaCP[timestep][ii]);
			}
		}
		if (inputScalars.storeResidual) {
			const float normi = static_cast<float>(af::norm(vec.pCP[timestep][subIter]));
			residual[kk] += (normi * normi * .5) + af::dot<float>(vec.pCP[timestep][subIter], y.as(f32));
		}
		vec.pCP[timestep][subIter] = input.copy();
		if (DEBUG) {
			mexPrintBase("w_vec.sigmaCP[timestep] = %f\n", w_vec.sigmaCP[timestep][ii]);
			mexEval();
		}
	}
	else if (MethodList.PDHGKL) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing PDHG with KL");
		if (w_vec.precondTypeMeas[0] || w_vec.precondTypeMeas[1]) {
			if (w_vec.precondTypeMeas[1]) {
				af::array apu1 = y.copy();
				status = applyMeasPreconditioning(w_vec, inputScalars, apu1, proj, timestep, subIter);
				if (status != 0)
					return -1;
				status = applyMeasPreconditioning(w_vec, inputScalars, input, proj, timestep, subIter);
				if (status != 0)
					return -1;
				apu1 = af::select(apu1 < 0.f, 0.f, apu1);
				input = af::select(input < 0.f, 0.f, input);
				input = .5f * (1.f + vec.pCP[timestep][subIter] + w_vec.sigmaCP[timestep][ii] * input - af::sqrt(af::pow(vec.pCP[timestep][subIter] + w_vec.sigmaCP[timestep][ii] * input - 1.f, 2.) + 4.f * w_vec.sigmaCP[timestep][ii] * apu1));
			}
			else {
				if (inputScalars.verbose >= 3)
					mexPrint("Applying diagonal normalization preconditioner (1 / (A1)), type 0");
				input = .5f * (1.f + vec.pCP[timestep][subIter] + w_vec.sigmaCP[timestep][ii] * input / w_vec.M[timestep][subIter] - af::sqrt(af::pow(vec.pCP[timestep][subIter] + w_vec.sigmaCP[timestep][ii] * input / w_vec.M[timestep][subIter] - 1.f, 2.) + 4.f * w_vec.sigmaCP[timestep][ii] * y.as(f32) / w_vec.M[timestep][subIter]));
			}
		}
		else
			input = .5f * (1.f + vec.pCP[timestep][subIter] + w_vec.sigmaCP[timestep][ii] * input - af::sqrt(af::pow(vec.pCP[timestep][subIter] + w_vec.sigmaCP[timestep][ii] * input - 1.f, 2.) + 4.f * w_vec.sigmaCP[timestep][ii] * y.as(f32)));
		vec.pCP[timestep][subIter] = input.copy();
	}
	else if (MethodList.PDHGL1) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing CPL1/TVL1/TGVL1");
		af::array res = input - y.as(f32);
		status = applyMeasPreconditioning(w_vec, inputScalars, res, proj, timestep, subIter);
		if (status != 0)
			return -1;
		input = (vec.pCP[timestep][subIter] + w_vec.sigmaCP[timestep][ii] * res);
		input /= (af::max)(1.f, af::abs(input));
		vec.pCP[timestep][subIter] = input.copy();
	}
	else if (MethodList.SAGA) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing SAGA");
		if (inputScalars.CT) {
			input = af::select(input < 0.f, 0.f, input);
			if (inputScalars.verbose >= 3)
				mexPrint("CT mode");
			if (inputScalars.verbose >= 3 && inputScalars.randoms_correction)
				mexPrint("Adding scatter data to forward projection");
			if (inputScalars.randoms_correction) {
				const af::array expInput = (hostGain ? af::array(aNorm * af::exp(-input)) : af::array(af::exp(-input)));
				input = expInput * inputScalars.flat - (y.as(f32) * expInput) / (expInput + randomsData);
			}
			else
				input = (hostGain ? af::array(aNorm * af::exp(-input)) : af::array(af::exp(-input))) * inputScalars.flat - y.as(f32);
		}
		else {
			if (inputScalars.verbose >= 3)
				mexPrint("PET/SPECT mode");
			// List-mode data: the -1 term is replaced by the sensitivity image D (subtracted in computeOSEstimates)
			if (inputScalars.listmode > 0)
				input = y.as(f32) / (input);
			else
				input = y.as(f32) / (input) - 1.f;
		}
		status = applyMeasPreconditioning(w_vec, inputScalars, input, proj, timestep, subIter);
		if (status != 0)
			return -1;
}
	if (MethodList.CPType && inputScalars.subsetsUsed > 1) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing PDHG with subsets");
		input -= vec.p0CP[timestep];
	}
	if (MethodList.FISTA || MethodList.FISTAL1) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing FISTA/L1");
		input -= y.as(f32);
		status = applyMeasPreconditioning(w_vec, inputScalars, input, proj, timestep, subIter);
		if (inputScalars.storeResidual) {
			residual[kk] = af::norm(input);
			residual[kk] = residual[kk] * residual[kk] * .5f;
		}
		if (status != 0)
			return -1;
	}
	input = af::select(af::isInf(input), inputScalars.epps, input);
	input.eval();
	if (DEBUG || inputScalars.verbose >= 3) {
		af::sync();
		proj.tEndLocal = std::chrono::steady_clock::now();
		const std::chrono::duration<double> tDiff = proj.tEndLocal - proj.tStartLocal;
		mexPrintBase("Computations required for backprojection completed in %f seconds\n", tDiff.count());
	}
	return 0;
}