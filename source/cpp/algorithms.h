#pragma once
#include "functions.hpp"

inline af::array EM(const af::array& im, const af::array& Summ, const af::array& rhs)
{
	return (im / Summ * rhs);
}

/// <summary>
/// This function computes one-step-late phase for COSEM/OSEM/MLEM-based reconstructions
/// </summary>
/// <param name="Summ Sensitivity image for the current subset"></param>
/// <param name="dU Gradient of the selected prior"></param>
/// <param name="beta Regularization parameter"></param>
/// <param name="epps Small value used to prevent division by zero"></param>
/// <returns>OSL regularized sensitivity image</returns>
inline af::array OSL(const af::array& Summ, const af::array& dU, const float epps)
{
	return (Summ + dU + epps);
}

inline int MBSREM(af::array& im, af::array& rhs, const float U, const float* lam, const uint32_t iter, const uint32_t osa_iter,
	const scalarStruct& inputScalars, Weighting& w_vec, ProjectorClass& proj, const int ii = 0)
{
	int status = 0;
	const uint32_t kk = iter * inputScalars.subsets + osa_iter;
	af::array output;
	const af::array pp = im >= (U / 2.f);
	//if (beta != 0.f)
	//	rhs -= beta * dU;
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
	//if (beta != 0.f)
	//	rhs -= beta * dU;
	if (DEBUG) {
		mexPrintBase("U = %f\n", U);
		mexPrintBase("iter = %d\n", iter);
		//mexPrintBase("beta = %f\n", beta);
		mexPrintBase("lam[iter] = %f\n", lam[iter]);
		mexPrintBase("w_vec.dP = %f\n", af::sum<float>(w_vec.dP[ii]));
		mexEval();
	}
	status = applyImagePreconditioning(w_vec, inputScalars, rhs, im, proj, kk, ii);
	if (status != 0)
		return -1;
	im += lam[iter] * rhs * w_vec.dP[ii];
	im(im <= 0) = inputScalars.epps;
	//im(im > U) = U;
	return 0;
}

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

inline af::array ROSEM(const af::array& im, const af::array& Summ, const af::array& rhs, const float* lam, const uint32_t iter)
{
	return (im + lam[iter] * im / Summ * (rhs - Summ));
}

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

inline af::array DRAMA(const af::array& im, const af::array& Summ, const af::array& rhs, const float* lam, const uint32_t iter, const uint32_t sub_iter, const uint32_t subsets)
{
	return (im + lam[iter * subsets + sub_iter] * im / Summ * rhs);
}

// MAP-phase for BSREM and ROSEM-MAP
inline af::array MAP(const af::array& im, const float lam, const af::array& dU, const float epps)
{
	//af::array output = im - beta * dU;
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

inline int PKMA(af::array& im, af::array& rhs, Weighting& w_vec, const scalarStruct& inputScalars,
	const uint32_t iter, const uint32_t osa_iter, ProjectorClass& proj, const int ii = 0) {
	int status = 0;
	const uint32_t kk = iter * inputScalars.subsets + osa_iter;
	//const af::array im_ = im.copy();
	//if (DEBUG) {
	//	mexPrintBase("im = %f\n", af::sum<float>(im));
	//	mexPrintBase("rhs = %f\n", af::sum<float>(rhs));
	//	//mexPrintBase("dU = %f\n", af::sum<float>(dU));
	//	mexEval();
	//}
	//if (beta != 0.f)
	//	rhs += beta * dU;
	applyImagePreconditioning(w_vec, inputScalars, rhs, im, proj, kk, ii);
	//if (ii > 0)
	//	status = proj.PoissonUpdate(im, rhs, inputScalars, w_vec.lambda[iter] * (55.f), inputScalars.epps, w_vec.alphaM[kk], ii);
	//else
	if (inputScalars.computeRelaxation) {
		if (kk == 0 && ii == 0) {
		//if (iter == 0 && ii == 0) {
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
		//if (kerroin < 1.5f && kerroin > 0.f)
		//	w_vec.lambda[iter] *= (kerroin / 1.5f);
		//if (iter == 0) {
			if (kerroin < 1.5f && kerroin > 0.f)
				w_vec.lambda[iter] *= (kerroin / 1.5f);
		//}
		//else {
		//	//if (kerroin2 < 2.f && kerroin > 0.f)
		//		w_vec.lambda[iter] *= (kerroin2 / 2.f);
		//}
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
	//const af::array im_ = im.copy();
	//af_print_mem_info("mem pre", -1);
	af::array im_apu = im - w_vec.lambda[iter] * rhs;
	//af_print_mem_info("mem imapu1", -1);
	////if (DEBUG) {
	////	mexPrintBase("im_apu = %f\n", af::sum<float>(im_apu));
	////	mexPrintBase("rhs2 = %f\n", af::sum<float>(rhs));
	////	mexEval();
	////}
	if (inputScalars.enforcePositivity)
		im_apu(im_apu < inputScalars.epps) = inputScalars.epps;
	//af_print_mem_info("mem imapu2", -1);
	////im_apu = (1.f - w_vec.alphaM[kk]) * im_ + w_vec.alphaM[kk] * (w_vec.sigma_PKMA[kk] * im_apu);
	im = (1.f - w_vec.alphaM[kk]) * im + w_vec.alphaM[kk] * im_apu;
	//if (DEBUG) {
	//	mexPrintBase("kerroin = %f\n", kerroin);
	//	//mexPrintBase("dU = %f\n", af::sum<float>(dU));
	//	mexEval();
	//}
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

inline void PDHG1(af::array& rhs, const scalarStruct& inputScalars, Weighting& w_vec, AF_im_vectors& vec, const uint32_t subIter = 0, const int ii = 0) {
	if (inputScalars.adaptiveType == 1)
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

inline int PDHG2(af::array& im, af::array& rhs, const scalarStruct& inputScalars, Weighting& w_vec, AF_im_vectors& vec, ProjectorClass& proj, const uint32_t iter, const uint32_t subIter = 0, const int ii = 0) {
	int status = 0;
	const uint32_t kk = iter * inputScalars.subsets + subIter;
	af::array im_old;
	if (inputScalars.adaptiveType == 1)
		im_old = im.copy();
	status = applyImagePreconditioning(w_vec, inputScalars, rhs, im, proj, kk, ii);
	if (status != 0)
		return -1;
	if (inputScalars.subsets > 1) {
		if (inputScalars.verbose >= 3)
			mexPrint("Using PDHG w/ subsets");
		//im -= w_vec.tauCP * rhs;
		//im.eval();
		//if (inputScalars.enforcePositivity)
		//	im(im < inputScalars.epps) = inputScalars.epps;
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
		//const af::array uPrev = vec.uCP.copy();
		//vec.uCP -= w_vec.tauCP * rhs;
		//vec.uCP.eval();
		//if (inputScalars.enforcePositivity)
		//	vec.uCP(vec.uCP < inputScalars.epps) = inputScalars.epps;
		//im = vec.uCP + w_vec.thetaCP[kk] * (vec.uCP - uPrev);
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
		//im = rhs.copy();
	}
	//if (DEBUG) {
	//	mexPrintBase("vec.im_os = %f\n", af::sum<float>(im));
	//	mexPrintBase("rhs = %f\n", af::sum<float>(rhs));
	//	mexPrintBase("w_vec.tauCP[ii] = %f\n", w_vec.tauCP[ii]);
	//	mexPrintBase("inputScalars.epps = %f\n", inputScalars.epps);
	//	mexEval();
	//}
	if (ii == 0 && inputScalars.adaptiveType == 1) {
		const af::array q = (im_old - im) / w_vec.tauCP[ii] + inputScalars.subsets * vec.rhsCP[ii];
		const float w = af::dot<float>((im_old - im), q) / (static_cast<float>(af::norm((im_old - im)) * af::norm(q)));
		//const float q = af::sum<float>(af::abs((im_old - im) / w_vec.tauCP[ii] + inputScalars.subsets * vec.rhsCP[ii]));
		//float w;
		//if (iter == 0 && subIter == 0)
		//	w = inputscalars.subsets * af::sum<float>(af::abs(-vec.pPrevCP / w_vec.sigmaCP[ii] + vec.fpCP[subIter]));
		//else
			//w = inputScalars.subsets * af::sum<float>(af::abs(-vec.pPrevCP / w_vec.sigmaCP[ii] - vec.fpCP[subIter]));
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
		if (inputScalars.verbose >= 3) {
			//mexPrintBase("w = %f\n", w);
			//mexPrintBase("q = %f\n", q);
			mexPrintBase("w_vec.alphaCP[ii] = %f\n", w_vec.alphaCP[ii]);
			mexPrintBase("w_vec.tauCP = %f\n", w_vec.tauCP[ii]);
			mexPrintBase("w_vec.sigmaCP = %f\n", w_vec.sigmaCP[ii]);
			mexEval();
		}
	}
	//if (q > w_vec.LCP[ii] * w * 1.01f) {
	//	w_vec.tauCP[ii] = w_vec.tauCP[ii] / (1.f - w_vec.alphaCP[ii]);
	//	w_vec.sigmaCP[ii] = w_vec.sigmaCP[ii] * (1.f - w_vec.alphaCP[ii]);
	//	w_vec.alphaCP[ii] *= 0.99f;
	//}
	//else if (q < w_vec.LCP[ii] * w / 1.01f) {
	//	w_vec.sigmaCP[ii] = w_vec.sigmaCP[ii] / (1.f - w_vec.alphaCP[ii]);
	//	w_vec.tauCP[ii] = w_vec.tauCP[ii] * (1.f - w_vec.alphaCP[ii]);
	//	w_vec.alphaCP[ii] *= 0.99f;
	//}
	//if (inputScalars.verbose >= 3) {
	//	//mexPrintBase("vec.uCP = %f\n", af::sum<float>(vec.uCP));
	//	//mexPrintBase("im = %f\n", af::sum<float>(im));
	//	//mexPrintBase("rhs = %f\n", af::sum<float>(rhs));
	//	mexPrintBase("w_vec.tauCP = %f\n", w_vec.tauCP[ii]);
	//	mexPrintBase("w_vec.LCP = %f\n", w_vec.LCP[ii]);
	//	mexPrintBase("w_vec.thetaCP[kk] = %f\n", w_vec.thetaCP[kk]);
	//	mexPrintBase("w_vec.sigmaCP = %f\n", w_vec.sigmaCP[ii]);
	//	mexEval();
	//}
	return status;
}

//int CPTV(af::array& im, af::array& rhs, const scalarStruct& inputScalars, Weighting& w_vec, AF_im_vectors& vec, ProjectorClass& proj, const uint32_t iter, const uint32_t subIter) {
//	int status = 0;
//	af::sync();
//	status = proj.CPTVHelperQ(vec.qProxTV, w_vec.betaReg, w_vec.UseL2Ball);
//	af::sync();
//	if (status != 0)
//		return -1;
//	status = proj.CPTVDiv(vec.qProxTV, rhs, inputScalars, w_vec.betaReg);
//	if (DEBUG) {
//		mexPrintBase("w_vec.betaReg = %f\n", w_vec.betaReg);
//		mexPrintBase("w_vec.sigma2CP = %f\n", w_vec.sigma2CP);
//		mexPrintBase("vec.qProxTV = %f\n", af::sum<float>(vec.qProxTV));
//		mexEval();
//	}
//	af::sync();
//	if (status != 0)
//		return -1;
//	CPLS(im, rhs, inputScalars, w_vec, vec, proj, iter, subIter);
//	status = proj.CPTVGrad(im, vec.qProxTV, inputScalars, w_vec.sigma2CP, vec.apu);
//	if (status != 0)
//		return -1;
//	if (DEBUG) {
//		mexPrintBase("vec.qCPTVgrad = %f\n", af::sum<float>(vec.qProxTV));
//		mexEval();
//	}
//	return 0;
//}
//
//int CPTGV(af::array& im, af::array& rhs, const scalarStruct& inputScalars, Weighting& w_vec, AF_im_vectors& vec, ProjectorClass& proj, const uint32_t iter, const uint32_t subIter) {
//	int status = 0;
//	status = proj.CPTVGrad(im, vec.qProxTV, inputScalars, w_vec.sigma2CP, vec.vProxTGV);
//	af::sync();
//	if (status != 0)
//		return -1;
//	status = proj.CPTVHelperQ(vec.qProxTV, w_vec.alpha0CPTGV, w_vec.UseL2Ball);
//	if (status != 0)
//		return -1;
//	if (DEBUG) {
//		mexPrintBase("vec.qProxTV = %f\n", af::sum<float>(vec.qProxTV));
//		mexEval();
//	}
//	status = proj.CPTGVSymmDeriv(vec.vProxTGV, vec.qProxTGV, vec.qCPTGV2, inputScalars, w_vec.sigma2CP);
//	af::sync();
//	if (status != 0)
//		return -1;
//	status = proj.CPTGVHelperQ(vec.qProxTGV, vec.qCPTGV2, w_vec.alpha1CPTGV, w_vec.UseL2Ball);
//	if (DEBUG) {
//		mexPrintBase("vec.qCPTGV2 = %f\n", af::sum<float>(vec.qCPTGV2));
//		mexEval();
//	}
//	af::sync();
//	if (status != 0)
//		return -1;
//	status = proj.CPTVDiv(vec.qProxTV, rhs, inputScalars, w_vec.alpha1CPTGV);
//	if (DEBUG) {
//		mexPrintBase("w_vec.alpha0CPTGV = %f\n", w_vec.alpha0CPTGV);
//		mexPrintBase("w_vec.alpha1CPTGV = %f\n", w_vec.alpha1CPTGV);
//		mexPrintBase("w_vec.sigma2CP = %f\n", w_vec.sigma2CP);
//		mexPrintBase("vec.qProxTGV = %f\n", af::sum<float>(vec.qProxTGV));
//		mexEval();
//	}
//	af::sync();
//	if (status != 0)
//		return -1;
//	CPLS(im, rhs, inputScalars, w_vec, vec, proj, iter, subIter);
//	status = proj.CPTGVDiv(vec.qProxTGV, vec.qCPTGV2, vec.vProxTGV, vec.qProxTV, inputScalars, w_vec.alpha0CPTGV, w_vec.thetaCP, w_vec.tauCP);
//	if (status != 0)
//		return -1;
//	if (DEBUG) {
//		mexPrintBase("vec.vProxTGV = %f\n", af::sum<float>(vec.vProxTGV));
//		mexEval();
//	}
//	return 0;
//}
//
//
//int CPGrad(af::array& im, af::array& rhs, const scalarStruct& inputScalars, Weighting& w_vec, AF_im_vectors& vec, ProjectorClass& proj, const uint32_t iter, const uint32_t subIter, 
//	const float beta, const af::array& dU) {
//	int status = 0;
//	rhs += beta * dU;
//	if (DEBUG) {
//		mexPrintBase("dU = %f\n", af::sum<float>(dU));
//		mexEval();
//	}
//	rhs.eval();
//	status = CPLS(im, rhs, inputScalars, w_vec, vec, proj, iter, subIter);
//	return status;
//}

inline int FISTA(af::array& im, af::array& rhs, const scalarStruct& inputScalars, Weighting& w_vec, AF_im_vectors& vec, ProjectorClass& proj, const uint32_t iter = 0, const int ii = 0) {
	int status = 0;
	af::array uPrev = im.copy();
	status = applyImagePreconditioning(w_vec, inputScalars, rhs, im, proj, iter, ii);
	if (status != 0)
		return -1;
	im = vec.uFISTA[ii] - w_vec.tauCP[ii] * rhs;
	if (ii == 0) {
		const uint32_t it = iter + 1;
		w_vec.betaFISTA = static_cast<float>(it - 1) / static_cast<float>(it + 2);
		if (w_vec.betaFISTA <= 0.f) {
			w_vec.tFISTA = (1.f + std::sqrt(1.f + 4.f * w_vec.tNFista * w_vec.tNFista)) / 2.f;
			w_vec.betaFISTA = (w_vec.tNFista - 1.f) / w_vec.tFISTA;
			w_vec.tNFista = w_vec.tFISTA;
		}
	}
	vec.uFISTA[ii] = im + w_vec.betaFISTA * (im - uPrev);
	vec.uFISTA[ii].eval();
	im.eval();
	rhs.eval();
	return 0;
}

inline int FISTAL1(af::array& im, af::array& rhs, const scalarStruct& inputScalars, Weighting& w_vec, AF_im_vectors& vec, const float beta, ProjectorClass& proj, const uint32_t iter = 0, const int ii = 0) {
	int status = 0;
	status = FISTA(im, rhs, inputScalars, w_vec, vec, proj, iter, ii);
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