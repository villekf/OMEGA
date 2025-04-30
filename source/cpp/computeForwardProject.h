#pragma once
#include "functions.hpp"

// Computes all computations using the forward projection and outputting a measurement-domain vector
inline int computeForwardStep(const RecMethods& MethodList, af::array& y, af::array& input, const int64_t length, const scalarStruct& inputScalars,
	Weighting& w_vec, const af::array& randomsData, AF_im_vectors& vec, ProjectorClass& proj, const uint32_t iter = 0, const uint32_t subIter = 0, const int ii = 0, 
	float* residual = nullptr) {

	af::deviceGC();
	int status = 0;
	af::array indeksit;
	bool indS = false;
	uint32_t kk = inputScalars.currentSubset + iter * inputScalars.subsets;
	if (DEBUG) {
		mexPrintBase("length = %d\n", length);
		mexPrintBase("y.dims(0) = %d\n", y.dims(0));
		mexPrintBase("input.dims(0) = %d\n", input.dims(0));
		mexEval();
	}
	if (!inputScalars.CT)
		input(input > 0.f) += inputScalars.epps;
	if (w_vec.precondTypeMeas[1] && w_vec.filterIter > 0 && kk == w_vec.filterIter) {
		if (inputScalars.verbose >= 3)
			mexPrint("Filter iterations complete. Switching tau/sigma-values");
		if (MethodList.CPType || MethodList.FISTA || MethodList.FISTAL1 || MethodList.ProxTGV || MethodList.ProxTV) {
			for (int ii = 0; ii <= inputScalars.nMultiVolumes; ii++) {
				if (w_vec.sigmaCP[ii] == 1.f)
					w_vec.tauCP[ii] = w_vec.tauCP2[ii];
				else if (w_vec.sigmaCP[ii] == w_vec.tauCP[ii]) {
					w_vec.tauCP[ii] = w_vec.tauCP2[ii];
					w_vec.sigmaCP[ii] = w_vec.tauCP2[ii];
				}
				else
					w_vec.sigmaCP[ii] = w_vec.tauCP2[ii];
				//if (MethodList.CPType)
				//	w_vec.sigma2CP[ii] = w_vec.sigmaCP[ii];
				if (inputScalars.adaptiveType == 1)
					w_vec.alphaCP[ii] = 1.f;
				else if (inputScalars.adaptiveType == 2)
					w_vec.alphaCP[ii] = .95f;
			}
			w_vec.LCP = w_vec.LCP2;
		}
		if (MethodList.MRAMLA || MethodList.MBSREM || MethodList.SPS || MethodList.RAMLA || MethodList.BSREM || MethodList.ROSEM || MethodList.ROSEMMAP || MethodList.PKMA || MethodList.SAGA)
			w_vec.lambda = w_vec.lambdaFiltered;
		w_vec.precondTypeMeas[1] = false;
	}
	if (inputScalars.randoms_correction) {
		if ((MethodList.MBSREM || MethodList.MRAMLA || MethodList.SPS) && !inputScalars.CT && !af::allTrue<bool>(randomsData > 0.f)) {
			indeksit = y > 0.f && randomsData == 0.f && input <= w_vec.epsilon_mramla;
			indS = af::anyTrue<bool>(indeksit);
		}
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
			if (inputScalars.TOF)
				for (int to = 0; to < inputScalars.nBins; to++)
					input(af::seq(randomsData.elements() * to, randomsData.elements() * (to + 1) - 1)) = input(af::seq(randomsData.elements() * to, randomsData.elements() * (to + 1) - 1)) + randomsData;
			else
				input += randomsData;
			input.eval();
		}
	}
	if (MethodList.CPType && inputScalars.subsetsUsed > 1) {
		vec.p0CP = vec.pCP[subIter].copy();
	}
	if (MethodList.ACOSEM || MethodList.OSLCOSEM > 0 || MethodList.OSEM || MethodList.COSEM || MethodList.ECOSEM ||
		MethodList.ROSEM || MethodList.OSLOSEM || MethodList.ROSEMMAP) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing (OSL-)(A)(E)(C)(R)OSEM");
		if (inputScalars.CT) {
			input(input < 0.f) = 0.f;
			if (inputScalars.verbose >= 3)
				mexPrint("CT mode");
			if (inputScalars.verbose >= 3 && inputScalars.randoms_correction)
				mexPrint("Adding scatter data to forward projection");
			if (MethodList.OSEM || MethodList.ROSEM || MethodList.OSLOSEM || MethodList.ROSEMMAP)
				input = af::exp(-input);
			else
				input = af::exp(-input) / y;
		}
		else {
			if (inputScalars.verbose >= 3)
				mexPrint("PET/SPECT mode");
			input = y / (input);
		}
		input.eval();
	}
	else if (MethodList.RAMLA || MethodList.BSREM || MethodList.RBI || MethodList.RBIOSL || MethodList.DRAMA) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing RAMLA/BSREM/RBI/DRAMA");
		if (inputScalars.CT) {
			input(input < 0.f) = 0.f;
			if (inputScalars.verbose >= 3)
				mexPrint("CT mode");
			if (inputScalars.verbose >= 3 && inputScalars.randoms_correction)
				mexPrint("Adding scatter data to forward projection");
			if (inputScalars.randoms_correction)
				input = af::exp(-input) * inputScalars.flat - (y * af::exp(-input)) / (af::exp(-input) + randomsData);
			else
				input = af::exp(-input) * inputScalars.flat - y;
		}
		else {
			if (inputScalars.verbose >= 3)
				mexPrint("PET/SPECT mode");
				input = y / (input) - 1.f;
		}
		input.eval();
	}
	else if (MethodList.PKMA) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing PKMA");
		if (inputScalars.CT) {
			input(input < 0.f) = 0.f;
			if (inputScalars.verbose >= 3)
				mexPrint("CT mode");
			if (inputScalars.verbose >= 3 && inputScalars.randoms_correction)
				mexPrint("Adding scatter data to forward projection");
			if (inputScalars.randoms_correction)
				input = (y * af::exp(-input)) / (af::exp(-input) + randomsData) - af::exp(-input) * inputScalars.flat;
			else
				input = y - af::exp(-input) * inputScalars.flat;
		}
		else {
			if (inputScalars.verbose >= 3)
				mexPrint("PET/SPECT mode");
				input = 1.f - y / (input);
		}
		input.eval();
		status = applyMeasPreconditioning(w_vec, inputScalars, input, proj, subIter);
		if (status != 0)
			return -1;
	}
	else if (MethodList.MBSREM || MethodList.MRAMLA || MethodList.SPS) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing MBSREM/MRAMLA/SPS");
		if (inputScalars.CT) {
			input(input < 0.f) = 0.f;
			if (inputScalars.verbose >= 3)
				mexPrint("CT mode");
			if (inputScalars.verbose >= 3 && inputScalars.randoms_correction)
				mexPrint("Adding scatter data to forward projection");
			if (inputScalars.randoms_correction)
				input = af::exp(-input) * inputScalars.flat - (y * af::exp(-input)) / (af::exp(-input) + randomsData);
			else
				input = af::exp(-input) * inputScalars.flat - y;
		}
		else {
			if (inputScalars.verbose >= 3)
				mexPrint("PET/SPECT mode");
			if (inputScalars.randoms_correction && indS) {
				input(indeksit) = y(indeksit) / w_vec.epsilon_mramla - 1.f - (y(indeksit) / (w_vec.epsilon_mramla * w_vec.epsilon_mramla)) * (input(indeksit) - w_vec.epsilon_mramla);
				input(!indeksit) = y(!indeksit) / (input(!indeksit) + inputScalars.epps) - 1.f;
			}
			else
				input = y / (input) - 1.f;
		}
		input.eval();
		status = applyMeasPreconditioning(w_vec, inputScalars, input, proj, subIter);
		if (status != 0)
			return -1;
	}
	else if (MethodList.LSQR) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing LSQR");
		input -= w_vec.alphaLSQR * y;
		w_vec.betaLSQR = af::norm(input);
		input = input / w_vec.betaLSQR;
		input.eval();
		y = input;
		y.eval();
	}
	else if (MethodList.CGLS) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing CGLS");
		const float normi = af::norm(input);
		w_vec.alphaCGLS = w_vec.gammaCGLS / (normi * normi);
		input = vec.rCGLS - w_vec.alphaCGLS * input;
		input.eval();
		vec.rCGLS = input;
		vec.rCGLS.eval();
	}else if (MethodList.BB) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing BB");

		input  = input -y;
	}
	else if (MethodList.SART || MethodList.POCS) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing SART or ASD-POCS");
		if (MethodList.POCS)
			vec.f0POCS = vec.im_os;
		input = y - input;
		if (inputScalars.storeResidual) {
			//residual[kk] = af::sum<float>(af::matmulTN(input, input)) * .5;
			residual[kk] = af::norm(input);
			residual[kk] = residual[kk] * residual[kk] * .5f;
		}
		input /= w_vec.M[subIter];
		input.eval();
	}
	else if (MethodList.PDHG || MethodList.PDDY) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing PDHG");
		if (DEBUG) {
			mexPrintBase("rdim0 = %u\n", y.dims(0));
			mexPrintBase("rdim1 = %u\n", y.dims(1));
			mexEval();
		}
		af::array res = input - y;
		if (DEBUG) {
			mexPrintBase("res = %f\n", af::sum<float>(res));
			mexPrintBase("max(res) = %f\n", af::max<float>(res));
			mexEval();
		}
		if (inputScalars.storeResidual) {
			//af::array ressa = res.copy();
			//status = applyMeasPreconditioning(w_vec, inputScalars, ressa, proj, subIter);
			//residual[kk] = af::sum<float>(af::matmulTN(res, ressa)) * .5;
			residual[kk] = af::norm(res);
			residual[kk] = residual[kk] * residual[kk] * .5f;
		}
		status = applyMeasPreconditioning(w_vec, inputScalars, res, proj, subIter);
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
				input = (vec.pCP[subIter] + w_vec.sigmaCP[ii] * res);
				if (inputScalars.subsetType == 5 || inputScalars.subsetType == 4) {
					if (inputScalars.subsetType == 4)
						input = af::moddims(input, inputScalars.nRowsD, input.elements() / inputScalars.nRowsD);
					else
						input = af::moddims(input, inputScalars.nColsD, input.elements() / inputScalars.nColsD);
				}
				else
					input = af::moddims(input, inputScalars.nRowsD, inputScalars.nColsD, input.elements() / (inputScalars.nRowsD * inputScalars.nColsD));
				if (inputScalars.adaptiveType >= 1 && ii == 0) {
					w_vec.Ffilter = af::ifft(w_vec.filter) * w_vec.sigmaCP[ii];
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
					mexPrintBase("vec.pCP[subIter].dim(0) = %d\n", vec.pCP[subIter].dims(0));
					mexPrintBase("res.dim(0) = %d\n", res.dims(0));
					mexPrintBase("res.elements() = %d\n", res.elements());
					mexPrintBase("apu.elements() = %d\n", apu.elements());
					mexPrintBase("res.elements() / w_vec.filter.elements() = %d\n", res.elements() / apu.elements());
					mexEval();
				}
				input = (vec.pCP[subIter] + w_vec.sigmaCP[ii] * res);
			}
			//if (inputScalars.storeResidual)
			//	residual[kk] += static_cast<float>(af::norm(vec.pCP[subIter]) * af::norm(vec.pCP[subIter]) * .5) + af::dot<float>(vec.pCP[subIter], y);
			input.eval();
		}
		else {
			if (MethodList.ProxTGV) {
				if (inputScalars.verbose >= 3)
					mexPrint("Computing Proximal TGV");
				input = (vec.pCP[subIter] + w_vec.sigmaCP[ii] * res) / (1.f + w_vec.sigmaCP[ii] * w_vec.betaReg);
			}
			else {
				if (inputScalars.verbose >= 3)
					mexPrint("Computing PDHG");
				input = (vec.pCP[subIter] + w_vec.sigmaCP[ii] * res) / (1.f + w_vec.sigmaCP[ii]);
			}
			input.eval();
		}
		vec.pCP[subIter] = input.copy();
		if (DEBUG) {
			mexPrintBase("w_vec.sigmaCP = %f\n", w_vec.sigmaCP[ii]);
			mexEval();
		}
	}
	else if (MethodList.PDHGKL) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing PDHG/CV with KL");
		if (w_vec.precondTypeMeas[0] || w_vec.precondTypeMeas[1]) {
			if (w_vec.precondTypeMeas[1]) {
				af::array apu1 = y.copy();
				status = applyMeasPreconditioning(w_vec, inputScalars, apu1, proj, subIter);
				if (status != 0)
					return -1;
				status = applyMeasPreconditioning(w_vec, inputScalars, input, proj, subIter);
				if (status != 0)
					return -1;
				apu1(apu1 < 0.f) = 0.f;
				input(input < 0.f) = 0.f;
				input = .5f * (1.f + vec.pCP[subIter] + w_vec.sigmaCP[ii] * input - af::sqrt(af::pow(vec.pCP[subIter] + w_vec.sigmaCP[ii] * input - 1.f, 2.) + 4.f * w_vec.sigmaCP[ii] * apu1));
			}
			else {
				if (inputScalars.verbose >= 3)
					mexPrint("Applying diagonal normalization preconditioner (1 / (A1)), type 0");
				input = .5f * (1.f + vec.pCP[subIter] + w_vec.sigmaCP[ii] * input / w_vec.M[subIter] - af::sqrt(af::pow(vec.pCP[subIter] + w_vec.sigmaCP[ii] * input / w_vec.M[subIter] - 1.f, 2.) + 4.f * w_vec.sigmaCP[ii] * y / w_vec.M[subIter]));
			}
		}
		else
			input = .5f * (1.f + vec.pCP[subIter] + w_vec.sigmaCP[ii] * input - af::sqrt(af::pow(vec.pCP[subIter] + w_vec.sigmaCP[ii] * input - 1.f, 2.) + 4.f * w_vec.sigmaCP[ii] * y));
		input.eval();
		vec.pCP[subIter] = input.copy();
	}
	else if (MethodList.PDHGL1) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing CPL1/TVL1/TGVL1");
		af::array res = input - y;
		status = applyMeasPreconditioning(w_vec, inputScalars, res, proj, subIter);
		if (status != 0)
			return -1;
		input = (vec.pCP[subIter] + w_vec.sigmaCP[ii] * res);
		input /= (af::max)(1.f, af::abs(input));
		input.eval();
		vec.pCP[subIter] = input.copy();
	}
	else if (MethodList.SAGA) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing SAGA");
		if (inputScalars.CT) {
			input(input < 0.f) = 0.f;
			if (inputScalars.verbose >= 3)
				mexPrint("CT mode");
			if (inputScalars.verbose >= 3 && inputScalars.randoms_correction)
				mexPrint("Adding scatter data to forward projection");
			if (inputScalars.randoms_correction)
				input = af::exp(-input) * inputScalars.flat - (y * af::exp(-input)) / (af::exp(-input) + randomsData);
			else
				input = af::exp(-input) * inputScalars.flat - y;
		}
		else {
			if (inputScalars.verbose >= 3)
				mexPrint("PET/SPECT mode");
			input = y / (input) - 1.f;
		}
		input.eval();
		status = applyMeasPreconditioning(w_vec, inputScalars, input, proj, subIter);
		if (status != 0)
			return -1;
}
	if (MethodList.CPType && inputScalars.subsetsUsed > 1) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing PDHG with subsets");
		input -= vec.p0CP;
		input.eval();
	}
	if (MethodList.FISTA || MethodList.FISTAL1) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing FISTA/L1");
		input -= y;
		status = applyMeasPreconditioning(w_vec, inputScalars, input, proj, subIter);
		if (inputScalars.storeResidual) {
			//residual[kk] = af::sum<float>(af::matmulTN(input, input)) * .5;
			residual[kk] = af::norm(input);
			residual[kk] = residual[kk] * residual[kk] * .5f;
		}
		if (status != 0)
			return -1;
	}
	input(af::isNaN(input)) = inputScalars.epps;
	input(af::isInf(input)) = inputScalars.epps;
	input.eval();

	if (inputScalars.verbose >= 3)
		mexPrint("Computations required for backprojection completed");
	af::sync();
	af::deviceGC();
	return 0;
}