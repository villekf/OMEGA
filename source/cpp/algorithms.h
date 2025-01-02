#pragma once
#include "functions.hpp"
#include "priors.h"

/// <summary>
/// Computes the final MLEM/OSEM update
/// </summary>
/// <param name="im Current image estimate"></param>
/// <param name="Summ Sensitivity image for the current subset"></param>
/// <param name="rhs Backprojection"></param>
/// <returns>MLEM/OSEM image</returns>
inline af::array EM(const af::array& im, const af::array& Summ, const af::array& rhs)
{
	return (im / Summ * rhs);
}

/// <summary>
/// This function computes one-step-late phase for COSEM/OSEM/MLEM-based reconstructions
/// </summary>
/// <param name="Summ Sensitivity image for the current subset"></param>
/// <param name="dU Gradient of the selected prior"></param>
/// <param name="epps Small value used to prevent division by zero"></param>
/// <returns>OSL regularized sensitivity image</returns>
inline af::array OSL(const af::array& Summ, const af::array& dU, const float epps)
{
	return (Summ + dU + epps);
}

// Computes the final MBSREM update
// Can use image-based preconditioners which are applied before the final update
inline int MBSREM(af::array& im, af::array& rhs, const float U, const float* lam, const uint32_t iter, const uint32_t osa_iter,
	const scalarStruct& inputScalars, Weighting& w_vec, ProjectorClass& proj, const int ii = 0)
{
	int status = 0;
	const uint32_t kk = iter * inputScalars.subsets + osa_iter;
	af::array output;
	const af::array pp = im >= (U / 2.f);
	if (af::anyTrue<bool>(pp)) {
		af::array apuIm = im;
		apuIm(pp) = U - apuIm(pp);
		applyImagePreconditioning(w_vec, inputScalars, rhs, apuIm, proj, kk);
	}
	else
		applyImagePreconditioning(w_vec, inputScalars, rhs, im, proj, kk, ii);
#ifndef CPU
	status = poissonUpdateAF(im, rhs, inputScalars, lam[iter], inputScalars.epps, U, proj, ii);
#else
	im = im + lam[iter] * rhs;
	im(im < inputScalars.epps) = inputScalars.epps;
	im(im >= U) = U - inputScalars.epps;
#endif
	return status;
}

// BSREM subset update
inline int BSREM(af::array& im, const af::array& rhs, const float* lam, const uint32_t iter, const scalarStruct& inputScalars, ProjectorClass& proj, const int ii = 0)
{
	int status = 0;
#ifndef CPU
	status = poissonUpdateAF(im, rhs, inputScalars, lam[iter], inputScalars.epps, 1.f, proj, ii);
#else
	 im = (im + lam[iter] * im * rhs);
#endif
	return status;
}

// Subset-based separable paraboidal surrogates
inline int SPS(af::array& im, af::array& rhs, const float U, const float* lam, const uint32_t iter, const uint32_t osa_iter,
	const scalarStruct& inputScalars, Weighting& w_vec, ProjectorClass& proj, const int ii = 0) {
	int status = 0;
	const uint32_t kk = iter * inputScalars.subsets + osa_iter;
	if (DEBUG) {
		mexPrintBase("U = %f\n", U);
		mexPrintBase("iter = %d\n", iter);
		mexPrintBase("lam[iter] = %f\n", lam[iter]);
		mexPrintBase("w_vec.dP = %f\n", af::sum<float>(w_vec.dP[ii]));
		mexEval();
	}
	status = applyImagePreconditioning(w_vec, inputScalars, rhs, im, proj, kk, ii);
	if (status != 0)
		return -1;
	im += lam[iter] * rhs * w_vec.dP[ii];
	im(im <= 0) = inputScalars.epps;
	return 0;
}

// ECOSEM update
inline af::array ECOSEM(const af::array& im, const af::array& D, const af::array& OSEM_apu, const af::array& COSEM_apu, const float epps)
{
	float alpha_eco = 1.f;
	af::array output = alpha_eco * OSEM_apu + (1.f - alpha_eco) * COSEM_apu;
	float eco_s1 = af::sum<float>(D * (-COSEM_apu * af::log(im + epps) + im));
	float eco_s2 = af::sum<float>(D * (-COSEM_apu * af::log(output + epps) + output));
	while (alpha_eco > 0.0096f && eco_s1 < eco_s2) {
		alpha_eco *= 0.9f;
		output = alpha_eco * OSEM_apu + (1.f - alpha_eco) * COSEM_apu;
		eco_s2 = af::sum<float>(D * (-COSEM_apu * af::log(output + epps) + output));
	}
	if (alpha_eco <= 0.0096f)
		output = COSEM_apu;
	return output;
}

// ROSEM update
inline af::array ROSEM(const af::array& im, const af::array& Summ, const af::array& rhs, const float* lam, const uint32_t iter)
{
	return (im + lam[iter] * im / Summ * (rhs - Summ));
}

// RBI update
inline af::array RBI(const af::array& im, const af::array& Summ, const af::array& rhs, const af::array& D = af::constant(0.f, 1, 1), const float beta = 0.f, const af::array& dU = af::constant(0.f, 1, 1))
{
	af::array output = im;
	if (beta == 0.f) {
		const float Summa = 1.f / af::max<float>(Summ / D);
		output += Summa * (im / D) * (rhs);
	}
	else {
		const float Summa = 1.f / af::max<float>((Summ + dU) / (D + dU));
		output += Summa * (im / (D + dU)) * (rhs - dU);
	}
	return output;
}

// DRAMA update
inline af::array DRAMA(const af::array& im, const af::array& Summ, const af::array& rhs, const float* lam, const uint32_t iter, const uint32_t sub_iter, const uint32_t subsets)
{
	return (im + lam[iter * subsets + sub_iter] * im / Summ * rhs);
}

// MAP-phase for BSREM and ROSEM-MAP
// Iteration update
inline af::array MAP(const af::array& im, const float lam, const af::array& dU, const float epps)
{
	af::array output = im - lam * im * dU;
	output(output < epps) = epps;
	return output;
}

inline af::array COSEM(const af::array& im, const af::array& C_co, const af::array& D, const float h, const uint32_t COSEM_TYPE)
{
	af::array output;
	if (COSEM_TYPE == 1) {
		output = af::pow(af::sum(C_co, 1) / D, h);
	}
	else {
		output = (af::sum(C_co, 1) / D);
	}
	return output;
}

// PKMA
// Can use image-based preconditioners which are applied before the final update
inline int PKMA(af::array& im, af::array& rhs, Weighting& w_vec, const scalarStruct& inputScalars,
	const uint32_t iter, const uint32_t osa_iter, ProjectorClass& proj, const int ii = 0) {
	int status = 0;
	const uint32_t kk = iter * inputScalars.subsets + osa_iter;
	applyImagePreconditioning(w_vec, inputScalars, rhs, im, proj, kk, ii);
	if (inputScalars.computeRelaxation) {
		if (kk == 0 && ii == 0) {
			w_vec.lambda[iter] = af::norm(im) / af::norm(rhs) * .25f;
			const float kerroin = af::norm(im) / af::norm(rhs * w_vec.lambda[iter]);
			const float kerroin2 = std::fabs(af::max<float>(im) / af::max<float>(rhs));
			const float kerroin3 = af::median<float>(im) / af::median<float>(rhs);
			if (DEBUG) {
				mexPrintBase("kerroin = %f\n", kerroin);
				mexPrintBase("kerroinMax = %f\n", kerroin2);
				mexPrintBase("kerroinMed = %f\n", kerroin3);
				mexEval();
			}
		}
		else if (iter > 0 && osa_iter == 0 && ii == 0)
			w_vec.lambda[iter] = w_vec.lambda[iter - 1] * (1.f / (iter / 35.f + 1.f));
		const float kerroin = af::norm(im) / af::norm(rhs * w_vec.lambda[iter]);
		const float kerroin2 = std::fabs(af::max<float>(im) / af::max<float>(rhs));
		const float kerroin3 = af::median<float>(im) / af::median<float>(rhs);
		const float kerroin4 = af::norm(im) / af::norm(rhs);
		const float kerroin5 = af::norm(im / rhs);
		const float kerroin6 = af::norm(im + rhs);
		const float kerroin7 = af::mean<float>(im) / af::mean<float>(rhs);
		if (DEBUG) {
			mexPrintBase("kerroin = %f\n", kerroin);
			mexPrintBase("kerroinMax = %f\n", kerroin2);
			mexPrintBase("kerroinMed = %f\n", kerroin3);
			mexPrintBase("kerroin4 = %f\n", kerroin4);
			mexPrintBase("kerroinJako = %f\n", kerroin5);
			mexPrintBase("kerroinSumma = %f\n", kerroin6);
			mexPrintBase("kerroinMean = %f\n", kerroin7);
			mexPrintBase("w_vec.lambda[iter] = %f\n", w_vec.lambda[iter]);
			mexEval();
		}
	}
	if (inputScalars.relaxScaling) {
		const float kerroin = af::norm(im) / af::norm(rhs * w_vec.lambda[iter]);
		const float kerroin2 = std::fabs(af::max<float>(im) / af::max<float>(rhs * w_vec.lambda[iter]));
		const float kerroin3 = af::median<float>(im) / af::median<float>(rhs * w_vec.lambda[iter]);
			if (kerroin < 1.5f && kerroin > 0.f)
				w_vec.lambda[iter] *= (kerroin / 1.5f);
		if (DEBUG) {
			mexPrintBase("kerroin = %f\n", kerroin);
			mexPrintBase("kerroinMax = %f\n", kerroin2);
			mexPrintBase("kerroinMed = %f\n", kerroin3);
			mexPrintBase("w_vec.lambda[iter] = %f\n", w_vec.lambda[iter]);
			mexEval();
		}
	}
#ifndef CPU
	status = poissonUpdateAF(im, rhs, inputScalars, w_vec.lambda[iter], inputScalars.epps, w_vec.alphaM[kk], proj, ii);
#else
	af::array im_apu = im - w_vec.lambda[iter] * rhs;
	if (inputScalars.enforcePositivity)
		im_apu(im_apu < inputScalars.epps) = inputScalars.epps;
	im = (1.f - w_vec.alphaM[kk]) * im + w_vec.alphaM[kk] * im_apu;
#endif
	return status;
}

inline void LSQR(const scalarStruct& inputScalars, Weighting& w_vec, const uint32_t iter, AF_im_vectors& vec, const int ii = 0) {
	if (iter == 0)
		vec.wLSQR[ii] = vec.im_os[ii];
	vec.im_os[ii] = vec.rhs_os[ii] - w_vec.betaLSQR * vec.im_os[ii];
	if (ii == inputScalars.nMultiVolumes) {
		af::array temp = vec.im_os[0];
		for (int ll = 1; ll <= inputScalars.nMultiVolumes; ll++)
			temp = af::join(0, vec.im_os[ll], temp);
		w_vec.alphaLSQR = af::norm(temp);
		for (int ll = 0; ll <= inputScalars.nMultiVolumes; ll++)
			vec.im_os[ll] = vec.im_os[ll] / w_vec.alphaLSQR;
		const float rho_ = sqrt(w_vec.rhoLSQR * w_vec.rhoLSQR + w_vec.betaLSQR * w_vec.betaLSQR);
		const float c = w_vec.rhoLSQR / rho_;
		const float s = w_vec.betaLSQR / rho_;
		w_vec.thetaLSQR = s * w_vec.alphaLSQR;
		w_vec.rhoLSQR = -c * w_vec.alphaLSQR;
		const float phi_ = c * w_vec.phiLSQR;
		w_vec.phiLSQR = s * w_vec.phiLSQR;
		for (int ll = 0; ll <= inputScalars.nMultiVolumes; ll++) {
			vec.fLSQR[ll] = (phi_ / rho_) * vec.wLSQR[ll] + vec.fLSQR[ll];
			vec.fLSQR[ll].eval();
			vec.wLSQR[ll] = vec.im_os[ll] - (w_vec.thetaLSQR / rho_) * vec.wLSQR[ll];
			vec.wLSQR[ll].eval();
			if (iter == inputScalars.Niter - 1)
				vec.im_os[ll] = vec.fLSQR[ll];
		}
	}
}

inline void CGLS(const scalarStruct& inputScalars, Weighting& w_vec, const uint32_t iter, AF_im_vectors& vec, const int ii = 0) {
	if (ii == inputScalars.nMultiVolumes) {
		float gamma_ = 0.f;
		for (int ll = 0; ll <= inputScalars.nMultiVolumes; ll++)
			gamma_ += af::sum<float>(vec.rhs_os[ll] * vec.rhs_os[ll]);
		const float beta = gamma_ / w_vec.gammaCGLS;
		for (int ll = 0; ll <= inputScalars.nMultiVolumes; ll++) {
			vec.fCGLS[ll] = vec.fCGLS[ll] + w_vec.alphaCGLS * vec.im_os[ll];
			vec.fCGLS[ll].eval();
			if (iter == inputScalars.Niter - 1)
				vec.im_os[ll] = vec.fCGLS[ll];
			else
				vec.im_os[ll] = vec.rhs_os[ll] + beta * vec.im_os[ll];
		}
		w_vec.gammaCGLS = gamma_;
	}
}

inline af::array SART(const af::array& im, const af::array& Summ, const af::array& rhs, const float lam) {
	return (im + lam * (rhs / Summ));
}

// Necessary PDHG step, when using subsets, before regularization is applied
inline void PDHG1(af::array& rhs, const scalarStruct& inputScalars, Weighting& w_vec, AF_im_vectors& vec, const uint32_t subIter = 0, const int ii = 0) {
	if (inputScalars.adaptiveType >= 1)
		vec.rhsCP[ii] = rhs.copy();
	if (inputScalars.subsets > 1) {
		if (DEBUG) {
			mexPrintBase("rhs = %f\n", af::sum<float>(rhs));
			mexPrintBase("vec.uCP[ii] = %f\n", af::sum<float>(vec.uCP[ii]));
			mexPrintBase("w_vec.thetaCP[subIter] = %f\n", w_vec.thetaCP[subIter]);
			mexEval();
		}
		if (inputScalars.verbose >= 3)
			mexPrint("Using PDHG w/ subsets");
		vec.uCP[ii] += rhs;
		vec.uCP[ii].eval();
		rhs = vec.uCP[ii] + (static_cast<float>(inputScalars.subsets) * w_vec.thetaCP[subIter]) * rhs;
	}
}

// Final PDHG update step
// Different methods for subset and non-subset cases
inline int PDHG2(af::array& im, af::array& rhs, scalarStruct& inputScalars, Weighting& w_vec, AF_im_vectors& vec, ProjectorClass& proj, const uint32_t iter, const uint32_t subIter = 0, const int ii = 0, 
	const int64_t* pituus = nullptr, const af::array& g = af::constant(0.f, 0), const uint64_t m_size = 1, const std::vector<int64_t>& length = std::vector<int64_t>(0)) {
	int status = 0;
	const uint32_t kk = iter * inputScalars.subsets + subIter;
	af::array im_old;
	if (inputScalars.adaptiveType >= 1)
		im_old = im.copy();
	if (ii == 0) {
		status = applyImagePreconditioning(w_vec, inputScalars, rhs, im, proj, kk, ii);
		if (status != 0)
			return -1;
	}
	if (inputScalars.subsets > 1) {
		if (inputScalars.verbose >= 3)
			mexPrint("Using PDHG w/ subsets");
#ifndef CPU
		status = PDHGUpdateAF(im, rhs, inputScalars, vec, inputScalars.epps, 1.f, w_vec.tauCP[ii], proj, ii);
#else
		im -= w_vec.tauCP[ii] * rhs;
		im.eval();
		if (inputScalars.enforcePositivity)
			im(im < inputScalars.epps) = inputScalars.epps;
#endif
	}
	else {
		if (inputScalars.verbose >= 3)
			mexPrint("Using PDHG W/O subsets");
#ifndef CPU
		status = PDHGUpdateAF(im, rhs, inputScalars, vec, inputScalars.epps, w_vec.thetaCP[kk], w_vec.tauCP[ii], proj, ii);
#else
		const af::array uPrev = vec.uCP[ii].copy();
		vec.uCP[ii] -= w_vec.tauCP[ii] * rhs;
		vec.uCP[ii].eval();
		if (inputScalars.enforcePositivity)
			vec.uCP[ii](vec.uCP[ii] < inputScalars.epps) = inputScalars.epps;
		im = vec.uCP[ii] + w_vec.thetaCP[kk] * (vec.uCP[ii] - uPrev);
#endif
	}
	if ((w_vec.precondTypeMeas[1] && subIter + inputScalars.subsets * iter >= w_vec.filterIter) || !w_vec.precondTypeMeas[1]) {
		if (ii == 0 && inputScalars.adaptiveType == 1) {
			const af::array q = (im_old - im) / w_vec.tauCP[ii] + inputScalars.subsets * vec.rhsCP[ii];
			const float w = af::dot<float>((im_old - im), q) / (static_cast<float>(af::norm((im_old - im)) * af::norm(q)));
			if (w < 0.f) {
				w_vec.tauCP[ii] = w_vec.tauCP[ii] / (1.f + w_vec.alphaCP[ii]);
				w_vec.sigmaCP[ii] = w_vec.sigmaCP[ii] * (1.f + w_vec.alphaCP[ii]);
				w_vec.alphaCP[ii] *= 0.99f;
			}
			else if (w >= .999f) {
				w_vec.sigmaCP[ii] = w_vec.sigmaCP[ii] / (1.f + w_vec.alphaCP[ii]);
				w_vec.tauCP[ii] = w_vec.tauCP[ii] * (1.f + w_vec.alphaCP[ii]);
				w_vec.alphaCP[ii] *= 0.99f;
			}
			w_vec.sigma2CP[ii] = w_vec.sigmaCP[ii];
			if (inputScalars.verbose >= 3) {
				mexPrintBase("w_vec.alphaCP[ii] = %f\n", w_vec.alphaCP[ii]);
				mexPrintBase("w_vec.tauCP = %f\n", w_vec.tauCP[ii]);
				mexPrintBase("w_vec.sigmaCP = %f\n", w_vec.sigmaCP[ii]);
				mexPrintBase("w = %f\n", w);
				mexEval();
			}
		}
		else if (ii == 0 && inputScalars.adaptiveType == 2) {
			const af::array apu = vec.im_os[ii].copy();
			vec.im_os[ii] = (im_old - im);
			const float q = af::sum<float>(af::abs((vec.im_os[ii]) / w_vec.tauCP[ii] + inputScalars.subsets * vec.rhsCP[ii]));
			af::array outputFP = af::constant(0.f, m_size * inputScalars.nBins);
			status = forwardProjectionAFOpenCL(vec, inputScalars, w_vec, outputFP, subIter, length, g, m_size, proj, ii, pituus);
			if (status != 0) {
				return status;
			}
			const float w = inputScalars.subsets * af::sum<float>(af::abs(-vec.adapTypeA / w_vec.sigmaCP[ii] - outputFP));
			if (q > w * 1.01f * std::sqrt(w_vec.LCP[ii])) {
				w_vec.tauCP[ii] = w_vec.tauCP[ii] / (1.f - w_vec.alphaCP[ii]);
				w_vec.sigmaCP[ii] = w_vec.sigmaCP[ii] * (1.f - w_vec.alphaCP[ii]);
				w_vec.alphaCP[ii] *= 0.99f;
			}
			else if (q < (w * std::sqrt(w_vec.LCP[ii])) / 1.01f) {
				w_vec.sigmaCP[ii] = w_vec.sigmaCP[ii] / (1.f - w_vec.alphaCP[ii]);
				w_vec.tauCP[ii] = w_vec.tauCP[ii] * (1.f - w_vec.alphaCP[ii]);
				w_vec.alphaCP[ii] *= 0.99f;
			}
			w_vec.sigma2CP[ii] = w_vec.sigmaCP[ii];
			vec.im_os[ii] = apu.copy();
		}
	}
	return status;
}

inline int FISTA(af::array& im, af::array& rhs, const scalarStruct& inputScalars, Weighting& w_vec, AF_im_vectors& vec, ProjectorClass& proj, const uint32_t iter = 0, const uint32_t osa_iter = 0, const int ii = 0) {
	int status = 0;
	const uint32_t kk = iter * inputScalars.subsets + osa_iter;
	status = applyImagePreconditioning(w_vec, inputScalars, rhs, im, proj, kk, ii);
	if (status != 0)
		return -1;
	if (inputScalars.subsets > 1 && osa_iter == inputScalars.subsets - 1) {
		im -= w_vec.tauCP[ii] * rhs;
		if (ii == 0) {
			if (inputScalars.FISTAType == 1) {
				w_vec.tFISTA = (1.f + std::sqrt(1.f + 4.f * w_vec.tNFista * w_vec.tNFista)) / 2.f;
				w_vec.betaFISTA = (1.f - w_vec.tNFista) / w_vec.tFISTA;
				w_vec.tNFista = w_vec.tFISTA;
			}
			else {
				const uint32_t it = iter + 1;
				w_vec.betaFISTA = static_cast<float>(it - 1) / static_cast<float>(it + 2);
				if (w_vec.betaFISTA <= 0.f) {
					w_vec.tFISTA = (1.f + std::sqrt(1.f + 4.f * w_vec.tNFista * w_vec.tNFista)) / 2.f;
					w_vec.betaFISTA = (w_vec.tNFista - 1.f) / w_vec.tFISTA;
					w_vec.tNFista = w_vec.tFISTA;
				}
			}
		}
		im.eval();
		vec.uFISTA[ii] = im + w_vec.betaFISTA * (im - vec.uFISTA[ii]);
		vec.uFISTA[ii].eval();
	}
	else if (inputScalars.subsets == 1) {
		af::array uPrev = im.copy();
		im = vec.uFISTA[ii] - w_vec.tauCP[ii] * rhs;
		if (ii == 0) {
			if (inputScalars.FISTAType == 1) {
				w_vec.tFISTA = (1.f + std::sqrt(1.f + 4.f * w_vec.tNFista * w_vec.tNFista)) / 2.f;
				w_vec.betaFISTA = (w_vec.tNFista - 1.f) / w_vec.tFISTA;
				w_vec.tNFista = w_vec.tFISTA;
			}
			else {
				const uint32_t it = iter + 1;
				w_vec.betaFISTA = static_cast<float>(it - 1) / static_cast<float>(it + 2);
				if (w_vec.betaFISTA <= 0.f) {
					w_vec.tFISTA = (1.f + std::sqrt(1.f + 4.f * w_vec.tNFista * w_vec.tNFista)) / 2.f;
					w_vec.betaFISTA = (w_vec.tNFista - 1.f) / w_vec.tFISTA;
					w_vec.tNFista = w_vec.tFISTA;
				}
			}
		}
		vec.uFISTA[ii] = im + w_vec.betaFISTA * (im - uPrev);
		vec.uFISTA[ii].eval();
	}
	else {
		im -= w_vec.tauCP[ii] * rhs;
	}
	im.eval();
	rhs.eval();
	return 0;
}


// FISTA with L1 regularization
inline int FISTAL1(af::array& im, af::array& rhs, const scalarStruct& inputScalars, Weighting& w_vec, AF_im_vectors& vec, const float beta, ProjectorClass& proj, const uint32_t iter = 0, const uint32_t osa_iter = 0, const int ii = 0) {
	int status = 0;
	status = FISTA(im, rhs, inputScalars, w_vec, vec, proj, iter, osa_iter, ii);
	if (status != 0)
		return -1;
	const float a = w_vec.tauCP[ii] * beta;
	if (DEBUG) {
		mexPrintBase("a = %f\n", a);
		mexEval();
	}
	im(af::abs(im) <= a) = 0.f;
	af::array apu = -af::sign(im);
	apu(apu == 0.f) = 1.f;
	im = apu * (af::max)(af::abs(im) - a, 0.f);
	return 0;
}

// ASD-POCS
inline void POCS(af::array& im, scalarStruct& inputScalars, Weighting& w_vec, const RecMethods& MethodList, AF_im_vectors& vec, ProjectorClass& proj, const af::array& mData, const af::array& g, 
	std::vector<int64_t>& length, const int64_t* pituus, const uint32_t osa_iter, const uint32_t iter, const int ii = 0) {
	im(im < 0.f) = 0.f;
	if (DEBUG)
		mexPrint("Computing ASD-POCS");
	bool subsets = true;
	if (inputScalars.subsets > 1)
		subsets = osa_iter < inputScalars.subsets - 1;
	if (iter < inputScalars.Niter - 1 && subsets) {
		uint64_t m_size = length[osa_iter];
		if ((inputScalars.CT || inputScalars.SPECT || inputScalars.PET) && inputScalars.listmode == 0)
			m_size = static_cast<uint64_t>(inputScalars.nRowsD) * static_cast<uint64_t>(inputScalars.nColsD) * length[osa_iter];
		af::array outputFP = af::constant(0.f, m_size * inputScalars.nBins);
		int status = 0;
		status = forwardProjectionAFOpenCL(vec, inputScalars, w_vec, outputFP, osa_iter, length, g, m_size, proj, ii, pituus);
		if (status != 0) {
			return;
		}
		const float dd = static_cast<float>(af::norm(outputFP - mData));
		const float dp = static_cast<float>(af::norm(im - vec.f0POCS[ii]));
		if (DEBUG) {
			mexPrintBase("dd = %f\n", dd);
			mexEval();
		}
		if (iter == 0 && osa_iter == 0)
			w_vec.dtvg = w_vec.alphaPOCS * dp;
		vec.f0POCS[ii] = im;
		if (DEBUG) {
			mexPrintBase("dp = %f\n", dp);
			mexEval();
		}
		if (ii == 0) {
			for (int kk = 0; kk < w_vec.ng; kk++) {
				status = applyPrior(vec, w_vec, MethodList, inputScalars, proj, w_vec.beta, osa_iter + inputScalars.subsets * iter);
				if (status != 0)
					return;
				vec.dU /= (af::norm(vec.dU) + inputScalars.epps);
				im -= w_vec.dtvg * vec.dU;
				af::eval(im);
				af::eval(vec.dU);
			}
			const float dg = static_cast<float>(af::norm(im - vec.f0POCS[ii]));
			if (dg > w_vec.rMaxPOCS && dd > w_vec.POCSepps)
				w_vec.dtvg *= w_vec.POCSalphaRed;
		}
	}
}