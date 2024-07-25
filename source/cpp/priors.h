#pragma once
#include "functions.hpp"

inline int MRP(const af::array& im, const uint32_t medx, const uint32_t medy, const uint32_t medz, const scalarStruct& inputScalars,
	ProjectorClass& proj, af::array& dU, const float beta, const bool med_no_norm = false) {
	int status = 0;
	af::array padd = padding(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], medx, medy, medz);
	//const af::dim4 dimmi(padd.dims(0), padd.dims(1), padd.dims(2));
	//af::array grad = af::constant(0.f, dimmi);
	af::array grad = af::constant(0.f, im.elements());
	//padd = af::flat(padd);
	//grad = af::flat(grad);
	status = MRPAF(padd, grad, inputScalars, proj, medx, medy, medz);
	if (status != 0)
		return -1;
	//grad = af::moddims(grad, dimmi);
	//grad = grad(af::seq(medx, inputScalars.Nx[0] + medx - 1), af::seq(medy, inputScalars.Ny[0] + medy - 1),
	//	af::seq(medz, inputScalars.Nz[0] + medz - 1));
	//grad = af::flat(grad) + inputScalars.epps;
	//if (DEBUG) {
	//	mexPrintBase("grad = %f\n", af::sum<float>(grad));
	//	mexPrintBase("im = %f\n", af::sum<float>(im));
	//	mexPrintBase("erotus = %f\n", af::sum<float>(im - grad));
	//}
	if (med_no_norm)
		dU += (im - grad) * beta;
	else
		dU += (im - (grad)) / (grad + inputScalars.epps) * beta;
	af::sync();
	if (DEBUG) {
		mexPrintBase("min(grad2) = %f\n", af::min<float>(grad));
		mexPrintBase("grad2 = %f\n", af::sum<float>(grad));
		mexPrintBase("min(dU) = %f\n", af::min<float>(dU));
		mexPrintBase("max(dU) = %f\n", af::max<float>(dU));
		mexPrintBase("min(im) = %f\n", af::min<float>(im));
		mexPrintBase("max(im) = %f\n", af::max<float>(im));
	}
	return status;
}

inline af::array Quadratic_prior(const af::array& im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const scalarStruct& inputScalars, const af::array& weights_quad)
{
	const af::array apu_pad = padding(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], Ndx, Ndy, Ndz);
	af::array grad;
	af::array weights = weights_quad;
	if (Ndz == 0 || inputScalars.Nz[0] == 1) {
		grad = af::convolve2(apu_pad, weights);
		grad = grad(af::seq(Ndx, inputScalars.Nx[0] + Ndx - 1), af::seq(Ndy, inputScalars.Ny[0] + Ndy - 1), af::span);
	}
	else {
		grad = af::convolve3(apu_pad, weights);
		grad = grad(af::seq(Ndx, inputScalars.Nx[0] + Ndx - 1), af::seq(Ndy, inputScalars.Ny[0] + Ndy - 1), af::seq(Ndz, inputScalars.Nz[0] + Ndz - 1));
	}
	grad = af::flat(grad);
	return grad;
}

inline af::array Huber_prior(const af::array& im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const scalarStruct& inputScalars, const af::array& weights_huber, const float delta)
{
	af::array grad = Quadratic_prior(im, Ndx, Ndy, Ndz, inputScalars, weights_huber);
	if (af::sum<dim_t>(delta >= af::abs(af::flat(grad))) == grad.elements() && af::sum<int>(af::flat(grad)) != 0 && inputScalars.verbose > 0)
		mexPrint("Delta value of Huber prior larger than all the pixel difference values\n");
	grad(grad > delta) = delta;
	grad(grad < -delta) = -delta;
	return grad;
}

inline af::array FMH(const af::array& im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const scalarStruct& inputScalars, const uint32_t inffi,
	const af::array& offsets, const af::array& fmh_weights, const uint32_t alku_fmh, const bool med_no_norm = false)
{
	af::array grad;
	af::array indeksi1;
	const af::array padd = af::flat(padding(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], Ndx, Ndy, Ndz));
	uint32_t luup;
	if (inputScalars.Nz[0] == 1 || Ndz == 0) {
		grad = af::constant(0.f, inputScalars.im_dim[0], 5);
		luup = 4;
	}
	else {
		grad = af::constant(0.f, inputScalars.im_dim[0], 14);
		luup = 13;
	}
	for (uint32_t ii = 0; ii < luup; ii++) {
		indeksi1 = af::flat(offsets(af::span, af::seq(Ndx * ii, offsets.dims(1) - Ndx * (ii)-1, alku_fmh / Ndx - ii)));
		af::array apu_pad = af::moddims(padd(indeksi1 + 0), inputScalars.im_dim[0], fmh_weights.dims(0));
		grad(af::span, ii) = af::matmul(apu_pad, fmh_weights(af::span, ii));
	}
	indeksi1 = offsets.col(alku_fmh);
	grad(af::span, af::end) = padd(indeksi1 + 0U);
	grad = af::median(grad, 1);
	if (med_no_norm)
		grad = im - grad;
	else
		grad = (im - grad) / (grad + inputScalars.epps);
	return grad;
}

inline af::array L_filter(const af::array& im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const scalarStruct& inputScalars,
	const af::array& offsets, const af::array& a_L, const bool med_no_norm = false)
{
	af::array grad;
	af::array apu_pad = af::flat(padding(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], Ndx, Ndy, Ndz));
	apu_pad = apu_pad(af::flat(offsets));
	apu_pad = af::sort(af::moddims(apu_pad, inputScalars.im_dim[0], a_L.dims(0)), 1);
	grad = af::sum(af::batchFunc(apu_pad, af::transpose(a_L), batchMul), 1);
	if (med_no_norm)
		grad = im - grad;
	else
		grad = (im - grad) / (grad + inputScalars.epps);
	return grad;
}

inline af::array Weighted_mean(const af::array& im, const uint32_t Ndx, const uint32_t Ndy, const uint32_t Ndz, const scalarStruct& inputScalars,
	const af::array& weighted_weights, const float w_sum, const uint32_t mean_type = 1, const bool med_no_norm = false)
{
	af::array grad = af::constant(0.f, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0]);
	const float wsum = af::sum<float>(af::flat(weighted_weights));
	if (mean_type == 1U) {
		af::array padd = padding(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], Ndx, Ndy, Ndz);
		if (Ndz == 0 || inputScalars.Nz[0] == 1)
			grad = af::convolve2(padd, weighted_weights / wsum);
		else
			grad = af::convolve3(padd, weighted_weights / wsum);
	}
	else if (mean_type == 2U) {
		af::array padd = padding(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], Ndx, Ndy, Ndz);
		if (Ndz == 0 || inputScalars.Nz[0] == 1)
			grad = 1.f / af::convolve2(1.f / padd, weighted_weights / wsum);
		else
			grad = 1.f / af::convolve3(1.f / padd, weighted_weights / wsum);
	}
	else if (mean_type == 3U) {
		af::array padd = padding(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], Ndx, Ndy, Ndz);
		if (Ndz == 0 || inputScalars.Nz[0] == 1)
			grad = af::exp(af::convolve2(af::log(padd), weighted_weights / wsum));
		else
			grad = af::exp(af::convolve3(af::log(padd), weighted_weights / wsum));
	}
	else if (mean_type == 4U) {
		grad = af::constant(0.f, im.dims(0));
		af::array padd = padding(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], Ndx * 2, Ndy * 2, Ndz * 2);
		af::array m = af::convolve3(padd, weighted_weights / wsum);
		if (Ndz == 0 || inputScalars.Nz[0] == 1) {
			m = m(af::seq(Ndx, inputScalars.Nx[0] + Ndx * 3 - 1), af::seq(Ndy, inputScalars.Ny[0] + Ndy * 3 - 1), af::span);
		}
		else {
			m = m(af::seq(Ndx, inputScalars.Nx[0] + Ndx * 3 - 1), af::seq(Ndy, inputScalars.Ny[0] + Ndy * 3 - 1), af::seq(Ndz, inputScalars.Nz[0] + Ndz * 3 - 1));
		}
		af::array mm = af::constant(0.f, im.dims(0), weighted_weights.elements());
		int jj = 0;
		for (int ll = 0; ll < weighted_weights.dims(2); ll++) {
			for (int kk = 0; kk < weighted_weights.dims(1); kk++) {
				for (int hh = 0; hh < weighted_weights.dims(0); hh++) {
					af::array apu = m(af::seq(hh, hh + inputScalars.Nx[0] - 1), af::seq(kk, kk + inputScalars.Ny[0] - 1), af::seq(ll, ll + inputScalars.Nz[0] - 1));
					mm(af::span, jj) = af::flat(apu);
					jj++;
				}
			}
		}
		grad = batchFunc(im, mm, batchDiv);
		grad(grad < 1e-6f) = 1e-6f;
		grad = batchFunc(af::log(grad), af::transpose(af::flat(weighted_weights)), batchMul);
		grad = af::sum(grad, 1);
		grad = af::flat(grad);
	}
	else if (mean_type == 5U) {
		af::array padd = padding(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], Ndx * 2, Ndy * 2, Ndz * 2);
		af::array m = 1.f / af::convolve3(1.f / padd, weighted_weights / wsum);
		if (Ndz == 0 || inputScalars.Nz[0] == 1) {
			m = m(af::seq(Ndx, inputScalars.Nx[0] + Ndx * 3 - 1), af::seq(Ndy, inputScalars.Ny[0] + Ndy * 3 - 1), af::span);
		}
		else {
			m = m(af::seq(Ndx, inputScalars.Nx[0] + Ndx * 3 - 1), af::seq(Ndy, inputScalars.Ny[0] + Ndy * 3 - 1), af::seq(Ndz, inputScalars.Nz[0] + Ndz * 3 - 1));
		}
		af::array mm = af::constant(0.f, im.dims(0), weighted_weights.elements());
		int jj = 0;
		for (int ll = 0; ll < weighted_weights.dims(2); ll++) {
			for (int kk = 0; kk < weighted_weights.dims(1); kk++) {
				for (int hh = 0; hh < weighted_weights.dims(0); hh++) {
					af::array apu = m(af::seq(hh, hh + inputScalars.Nx[0] - 1), af::seq(kk, kk + inputScalars.Ny[0] - 1), af::seq(ll, ll + inputScalars.Nz[0] - 1));
					mm(af::span, jj) = af::flat(apu);
					jj++;
				}
			}
		}
		af::array im_t = im;
		im_t(im_t < 1e-6f) = 1e-6f;
		grad = batchFunc(batchFunc(im_t, mm, batchMinus), im_t, batchDiv);
		grad = grad - af::pow(batchFunc((batchFunc(im_t, mm, batchMinus)), std::sqrt(2.f) * (im_t), batchDiv), 2.);
		grad = batchFunc(grad, af::transpose(af::flat(weighted_weights)), batchMul);
		grad = af::sum(grad, 1);
		grad = af::flat(grad);
	}
	else if (mean_type == 6U) {
		af::array padd = padding(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], Ndx * 2, Ndy * 2, Ndz * 2);
		af::array m = af::exp(af::convolve3(af::log(padd), weighted_weights / wsum));
		if (Ndz == 0 || inputScalars.Nz[0] == 1) {
			m = m(af::seq(Ndx, inputScalars.Nx[0] + Ndx * 3 - 1), af::seq(Ndy, inputScalars.Ny[0] + Ndy * 3 - 1), af::span);
		}
		else {
			m = m(af::seq(Ndx, inputScalars.Nx[0] + Ndx * 3 - 1), af::seq(Ndy, inputScalars.Ny[0] + Ndy * 3 - 1), af::seq(Ndz, inputScalars.Nz[0] + Ndz * 3 - 1));
		}
		af::array mm = af::constant(0.f, im.dims(0), weighted_weights.elements());
		int jj = 0;
		for (int ll = 0; ll < weighted_weights.dims(2); ll++) {
			for (int kk = 0; kk < weighted_weights.dims(1); kk++) {
				for (int hh = 0; hh < weighted_weights.dims(0); hh++) {
					af::array apu = m(af::seq(hh, hh + inputScalars.Nx[0] - 1), af::seq(kk, kk + inputScalars.Ny[0] - 1), af::seq(ll, ll + inputScalars.Nz[0] - 1));
					mm(af::span, jj) = af::flat(apu);
					jj++;
				}
			}
		}
		af::array im_t = im;
		im_t(im_t < 1e-6f) = 1e-6f;
		grad = 1.f - batchFunc(mm, im_t, batchDiv);
		grad = batchFunc(grad, af::transpose(af::flat(weighted_weights)), batchMul);
		grad = af::sum(grad, 1);
		grad = af::flat(grad);
	}
	else {
		mexWarning("Unsupported mean type");
	}
	if (mean_type <= 3U) {
		if (Ndz == 0 || inputScalars.Nz[0] == 1) {
			grad = grad(af::seq(Ndx, inputScalars.Nx[0] + Ndx - 1), af::seq(Ndy, inputScalars.Ny[0] + Ndy - 1), af::span);
		}
		else {
			grad = grad(af::seq(Ndx, inputScalars.Nx[0] + Ndx - 1), af::seq(Ndy, inputScalars.Ny[0] + Ndy - 1), af::seq(Ndz, inputScalars.Nz[0] + Ndz - 1));
		}
		grad = af::flat(grad);
		if (med_no_norm)
			grad = im - grad;
		else
			grad = (im - grad) / (grad + inputScalars.epps);
	}
	return grad;
}

inline af::array AD(const af::array& im, const scalarStruct& inputScalars, const float TimeStepAD, const float KAD, const uint32_t NiterAD,
	const af_flux_function FluxType = AF_FLUX_EXPONENTIAL, const af_diffusion_eq DiffusionType = AF_DIFFUSION_GRAD, const bool med_no_norm = false)
{
	const af::array padd = af::moddims(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0]);
	af::array grad = af::anisotropicDiffusion(padd, TimeStepAD, KAD, NiterAD, FluxType, DiffusionType);
	grad = af::flat(grad);
	if (med_no_norm)
		grad = im - grad;
	else
		grad = (im - grad) / (grad + inputScalars.epps);
	return grad;
}


// Compute the TV prior
inline int TVprior(const scalarStruct& inputScalars, const TVdata& S, const af::array& ima, const Weighting& w_vec, ProjectorClass& proj, af::array& dU, 
	const float beta) {
	af::array gradi;
	int status = 0;

	status = TVAF(dU, ima, inputScalars, w_vec.data.SATVPhi, S, proj, beta);

	return status;
}

// Compute the hyperbolic prior
inline int hyperbolic(const scalarStruct& inputScalars, const af::array& ima, const Weighting& w_vec, ProjectorClass& proj, af::array& dU, const float beta) {
	af::array gradi;
	int status = 0;

	status = hyperAF(dU, ima, inputScalars, w_vec.data.SATVPhi, proj, beta);

	return status;
}

inline int proxTV(const af::array& im, const scalarStruct& inputScalars, AF_im_vectors& vec, ProjectorClass& proj, const Weighting& w_vec, af::array& dU, const float beta) {
	int status = 0;
#ifndef CPU
	status = proxTVGradAF(im, vec.qProxTV, inputScalars, w_vec.sigma2CP[0], vec.vProxTGV, proj);
	af::sync();
	if (status != 0)
		return -1;
	status = proxTVQAF(vec.qProxTV, beta, proj);
	af::sync();
	if (status != 0)
		return -1;
	status = proxTVDivAF(vec.qProxTV, dU, inputScalars, proj);
	if (DEBUG) {
		mexPrintBase("beta = %f\n", beta);
		mexPrintBase("w_vec.sigma2CP = %f\n", w_vec.sigma2CP[0]);
		mexPrintBase("vec.qProxTV = %f\n", af::sum<float>(vec.qProxTV[0]));
		mexEval();
	}
#else
	mexPrint("Proximal TV not supported with CPU implementation!");
	status = -1;
#endif
	return status;
}

inline int proxTGV(const af::array& im, const scalarStruct& inputScalars, AF_im_vectors& vec, ProjectorClass& proj, const Weighting& w_vec, af::array& dU, const uint32_t osa_iter = 0) {
	int status = 0;
#ifndef CPU
	status = proxTV(im, inputScalars, vec, proj, w_vec, dU, w_vec.alpha0CPTGV);
	if (DEBUG) {
		mexPrintBase("vec.qProxTV = %f\n", af::sum<float>(vec.qProxTV[0]));
		mexPrintBase("vec.qProxTGV = %f\n", af::sum<float>(vec.qProxTGV[0]));
		mexEval();
	}
	status = proxTGVSymmDerivAF(vec.vProxTGV, vec.qProxTGV, inputScalars, w_vec.sigma2CP[0], proj);
	af::sync();
	if (status != 0)
		return -1;
	status = proxTGVQAF(vec.qProxTGV, inputScalars, w_vec.alpha1CPTGV, proj);
	if (DEBUG) {
		mexPrintBase("vec.qCPTGV2 = %f\n", af::sum<float>(vec.qProxTGV[2]));
		mexEval();
	}
	af::sync();
	if (status != 0)
		return -1;
	//status = proj.ProxTVDiv(vec.qProxTV, dU, inputScalars);
	//af::sync();
	//if (status != 0)
	//	return -1;
	if (DEBUG) {
		mexPrintBase("w_vec.alpha0CPTGV = %f\n", w_vec.alpha0CPTGV);
		mexPrintBase("w_vec.alpha1CPTGV = %f\n", w_vec.alpha1CPTGV);
		mexPrintBase("w_vec.sigma2CP = %f\n", w_vec.sigma2CP);
		mexPrintBase("osa_iter = %d\n", osa_iter);
		mexPrintBase("vec.qProxTGV0 = %f\n", af::sum<float>(vec.qProxTGV[0]));
		mexEval();
	}
	status = proxTGVDivAF(vec.qProxTGV, vec.vProxTGV, vec.qProxTV, inputScalars, w_vec.thetaCP[osa_iter], w_vec.tauCP[0], proj);
#else
	mexPrint("Proximal TGV not supported with CPU implementation!");
	status = -1;
#endif
	return status;
}
//
//inline int proxRDP(const af::array& im, const scalarStruct& inputScalars, AF_im_vectors& vec, ProjectorClass& proj, af::array& dU, const float beta, const float gamma) {
//	int status = 0;
//	if (DEBUG) {
//		mexPrintBase("vec.qProx = %f\n", af::sum<float>(vec.qProx[0]));
//		mexPrintBase("im = %f\n", af::sum<float>(im));
//		mexEval();
//	}
//	status = proxRDP(im, vec.qProx, inputScalars, gamma);
//	af::sync();
//	if (status != 0)
//		return -1;
//	if (DEBUG) {
//		mexPrintBase("vec.qProx = %f\n", af::sum<float>(vec.qProx[0]));
//		mexEval();
//	}
//	status = proj.ProxHelperQ(vec.qProx, beta);
//	af::sync();
//	if (status != 0)
//		return -1;
//	if (DEBUG) {
//		mexPrintBase("w_vec.betaReg = %f\n", beta);
//		mexPrintBase("w_vec.gamma = %f\n", gamma);
//		mexPrintBase("vec.qProx = %f\n", af::sum<float>(vec.qProx[0]));
//		mexEval();
//	}
//	status = proj.ProxRDPTrans(vec.qProx, dU, inputScalars, gamma);
//	af::sync();
//	if (status != 0)
//		return -1;
//	return status;
//}
//
//inline int proxNLM(const af::array& im, const scalarStruct& inputScalars, AF_im_vectors& vec, ProjectorClass& proj, af::array& dU, const Weighting& w_vec, const float beta) {
//	int status = 0;
//	status = proj.ProxNLM(im, vec.qProx, inputScalars, w_vec);
//	af::sync();
//	if (status != 0)
//		return -1;
//	if (DEBUG) {
//		mexPrintBase("vec.qProx = %f\n", af::sum<float>(vec.qProx[0]));
//		mexEval();
//	}
//	status = proj.ProxHelperQ(vec.qProx, beta);
//	af::sync();
//	if (status != 0)
//		return -1;
//	if (DEBUG) {
//		mexPrintBase("w_vec.betaReg = %f\n", beta);
//		mexPrintBase("vec.qProx = %f\n", af::sum<float>(vec.qProx[0]));
//		mexEval();
//	}
//	status = proj.ProxRDPTrans(vec.qProx, dU, inputScalars, 0.f);
//	af::sync();
//	if (status != 0)
//		return -1;
//	return status;
//}

//af::array TGV(const af::array& im, const scalarStruct& inputScalars, const uint32_t maxits, const float alpha, const float beta)
//{
//	af::array grad = af::constant(0.f, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], f32);
//	af::array u = grad;
//	grad = af::flat(grad);
//	const af::array imi = af::moddims(im, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0]);
//	const float sigma = 1.f / 16.f;
//	const float tau = 1.f / 8.f;
//	af::array p1 = af::constant(0.f, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], f32);
//	af::array p2 = af::constant(0.f, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], f32);
//	af::array p3 = af::constant(0.f, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], f32);
//	af::array v1 = af::constant(0.f, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], f32);
//	af::array v2 = af::constant(0.f, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], f32);
//	af::array v3 = af::constant(0.f, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], f32);
//	af::array q1 = af::constant(0.f, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], f32);
//	af::array q2 = af::constant(0.f, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], f32);
//	af::array q3 = af::constant(0.f, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], f32);
//	af::array q4 = af::constant(0.f, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0], f32);
//
//	af::array vb1 = v1;
//	af::array vb2 = v2;
//	af::array vb3 = v3;
//	v1 = af::flat(v1);
//	v2 = af::flat(v2);
//	v3 = af::flat(v3);
//
//	for (uint32_t kk = 0; kk < maxits; kk++) {
//		af::array temp_u = u + imi;
//		af::array eta1 = af::join(0, af::diff1(temp_u), temp_u(0, af::span, af::span) - temp_u(af::end, af::span, af::span)) - vb1;
//		af::array eta2 = af::join(1, af::diff1(temp_u, 1), temp_u(af::span, 0, af::span) - temp_u(af::span, af::end, af::span)) - vb2;
//		af::array eta3 = af::join(2, af::diff1(temp_u, 2), temp_u(af::span, af::span, 0) - temp_u(af::span, af::span, af::end)) - vb3;
//		eta1 = p1 + sigma * eta1;
//		eta2 = p2 + sigma * eta2;
//		eta3 = p3 + sigma * eta3;
//		af::array apu = (af::max)(1.f, af::sqrt(eta1 * eta1 + eta2 * eta2 + eta3 * eta3) / beta);
//		apu = af::flat(apu);
//		p1 = af::flat(eta1) / (apu);
//		p1 = af::moddims(p1, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0]);
//		p2 = af::flat(eta2) / (apu);
//		p2 = af::moddims(p2, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0]);
//		p3 = af::flat(eta3) / (apu);
//		p3 = af::moddims(p3, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0]);
//		eta1 = af::join(0, af::diff1(vb1), vb1(0, af::span, af::span) - vb1(af::end, af::span, af::span));
//		eta2 = af::join(1, af::diff1(vb2, 1), vb2(af::span, 0, af::span) - vb2(af::span, af::end, af::span));
//		eta3 = af::join(2, af::diff1(vb3, 2), vb3(af::span, af::span, 0) - vb3(af::span, af::span, af::end));
//		af::array eta4 = (af::join(0, af::diff1(vb2), vb2(0, af::span, af::span) - vb2(af::end, af::span, af::span)) +
//			af::join(1, af::diff1(vb1, 1), vb1(af::span, 0, af::span) - vb1(af::span, af::end, af::span)) +
//			af::join(2, af::diff1(vb2, 2), vb2(af::span, af::span, 0) - vb2(af::span, af::span, af::end)) +
//			af::join(2, af::diff1(vb1, 2), vb1(af::span, af::span, 0) - vb1(af::span, af::span, af::end)) +
//			af::join(0, af::diff1(vb3), vb3(0, af::span, af::span) - vb3(af::end, af::span, af::span)) +
//			af::join(1, af::diff1(vb3, 1), vb3(af::span, 0, af::span) - vb3(af::span, af::end, af::span))) / 6.f;
//		eta1 = q1 + sigma * eta1;
//		eta2 = q2 + sigma * eta2;
//		eta3 = q3 + sigma * eta3;
//		eta4 = q4 + sigma * eta4;
//		apu = (af::max)(1.f, af::sqrt(eta1 * eta1 + eta2 * eta2 + eta3 * eta3 + eta4 * eta4) / alpha);
//		apu = af::flat(apu);
//		q1 = af::flat(eta1) / (apu);
//		q1 = af::moddims(q1, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0]);
//		q2 = af::flat(eta2) / (apu);
//		q2 = af::moddims(q2, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0]);
//		q3 = af::flat(eta3) / (apu);
//		q3 = af::moddims(q3, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0]);
//		q4 = af::flat(eta4) / (apu);
//		q4 = af::moddims(q4, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0]);
//
//		eta4 = af::join(0, p1(af::end, af::span, af::span) - p1(0, af::span, af::span), -af::diff1(p1)) +
//			af::join(1, p2(af::span, af::end, af::span) - p2(af::span, 0, af::span), -af::diff1(p2, 1)) +
//			af::join(2, p3(af::span, af::span, af::end) - p3(af::span, af::span, 0), -af::diff1(p3, 2));
//		af::array uold = grad;
//		grad -= tau * af::flat(eta4);
//		eta1 = af::join(0, q1(af::end, af::span, af::span) - q1(0, af::span, af::span), -af::diff1(q1)) +
//			af::join(1, q4(af::span, af::end, af::span) - q4(af::span, 0, af::span), -af::diff1(q4, 1)) +
//			af::join(2, q4(af::span, af::span, af::end) - q4(af::span, af::span, 0), -af::diff1(q4, 2)) - p1;
//		eta2 = af::join(0, q4(af::end, af::span, af::span) - q4(0, af::span, af::span), -af::diff1(q4)) +
//			af::join(1, q2(af::span, af::end, af::span) - q2(af::span, 0, af::span), -af::diff1(q2, 1)) +
//			af::join(2, q4(af::span, af::span, af::end) - q4(af::span, af::span, 0), -af::diff1(q4, 2)) - p2;
//		eta3 = af::join(0, q4(af::end, af::span, af::span) - q4(0, af::span, af::span), -af::diff1(q4)) +
//			af::join(2, q3(af::span, af::span, af::end) - q3(af::span, af::span, 0), -af::diff1(q3, 2)) +
//			af::join(1, q4(af::span, af::end, af::span) - q4(af::span, 0, af::span), -af::diff1(q4, 1)) - p3;
//		af::array v1old = v1;
//		af::array v2old = v2;
//		af::array v3old = v3;
//		v1 -= tau * af::flat(eta1);
//		v2 -= tau * af::flat(eta2);
//		v3 -= tau * af::flat(eta3);
//		u = 2.f * grad - uold;
//		u = af::moddims(u, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0]);
//		vb1 = af::moddims(2.f * v1 - v1old, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0]);
//		vb2 = af::moddims(2.f * v2 - v2old, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0]);
//		vb3 = af::moddims(2.f * v3 - v3old, inputScalars.Nx[0], inputScalars.Ny[0], inputScalars.Nz[0]);
//
//	}
//
//	return -grad;
//}

inline int RDP(const af::array& im, const scalarStruct& inputScalars, const float gamma, ProjectorClass& proj, af::array& dU, const float beta, const af::array& RDPref, 
	const bool RDPLargeNeighbor = false, const bool useRDPRef = false) {

	//af::array grad = af::constant(0.f, im.elements(), 1);
	int status = 0;
	af::sync();
	if (DEBUG) {
		mexPrintBase("im_RDP = %f\n", af::sum<float>(im));
		mexPrintBase("isnan(im_RDP) = %d\n", af::anyTrue<bool>(af::isNaN(im)));
		mexEval();
	}
	status = RDPAF(dU, im, inputScalars, gamma, proj, beta, RDPref, RDPLargeNeighbor, useRDPRef);
	if (DEBUG) {
		mexPrintBase("grad = %f\n", af::sum<float>(dU));
		mexPrintBase("min(grad) = %f\n", af::min<float>(dU));
		mexPrintBase("max(grad) = %f\n", af::max<float>(dU));
		mexEval();
	}
	return status;
}

inline int GGMRF(const af::array& im, const scalarStruct& inputScalars, const float p, const float q, const float c, const float pqc, ProjectorClass& proj, af::array& dU, const float beta) {

	//af::array grad = af::constant(0.f, im.elements(), 1);
	int status = 0;
	af::sync();
	if (DEBUG) {
		mexPrintBase("im_GGMRF = %f\n", af::sum<float>(im));
		mexPrintBase("isnan(im_GGMRF) = %d\n", af::anyTrue<bool>(af::isNaN(im)));
		mexEval();
	}
	status = GGMRFAF(dU, im, inputScalars, p, q, c, pqc, proj, beta);
	if (DEBUG) {
		mexPrintBase("grad = %f\n", af::sum<float>(dU));
		mexEval();
	}
	return status;
}



inline int NLM(ProjectorClass& proj, const af::array& im, Weighting& w_vec, const scalarStruct& inputScalars, af::array& dU, const float beta)
{
	int status = 0;
	int32_t type = 0;
	if (w_vec.NLTV)
		type = 1;
	else if (w_vec.NLM_MRP)
		type = 2;
	else
		type = 0;
	//af::array grad = af::constant(0.f, im.elements(), 1);
	af::sync();
	status = NLMAF(dU, im, inputScalars, w_vec, proj, beta);
	return status;
}


inline int applyPrior(AF_im_vectors& vec, Weighting& w_vec, const RecMethods& MethodList, const scalarStruct& inputScalars, ProjectorClass& proj, const float beta, const uint32_t osa_iter = 0, const uint8_t compute_norm_matrix = 0, const bool iter = false) {
	af::array* dU = nullptr;
	int status = 0;
	//if (MethodList.OSLOSEM) {
	//	if (compute_norm_matrix == 1u) {
	//		dU = &vec.Summ[0][0];
	//	}
	//	else if (compute_norm_matrix == 2u) {
	//		dU = &vec.Summ[0][osa_iter];
	//	}
	//	else
	//		return -1;
	//}
	//else if (MethodList.OSLCOSEM)
	//if (MethodList.OSLCOSEM)
	//	dU = &w_vec.D[0];
	if (iter)
		dU = &vec.im_os[0];
	else if (MethodList.RBIOSL || MethodList.OSLOSEM || MethodList.OSLCOSEM || MethodList.POCS) {
		vec.dU = af::constant(0.f, vec.im_os[0].elements());
		dU = &vec.dU;
	}
	else
		dU = &vec.rhs_os[0];
	if (MethodList.MRP) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing MRP gradient");
		status = MRP(vec.im_os[0], w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, proj, *dU, beta, w_vec.med_no_norm);
	}
	else if (MethodList.Quad) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing quadratic prior gradient");
		*dU += beta * Quadratic_prior(vec.im_os[0], w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.weights_quad);
	}
	else if (MethodList.Huber) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing Huber prior gradient");
		*dU += beta * Huber_prior(vec.im_os[0], w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.weights_huber, w_vec.huber_delta);
	}
	else if (MethodList.L) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing L-filter gradient");
		*dU += beta * L_filter(vec.im_os[0], w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.tr_offsets,
			w_vec.a_L, w_vec.med_no_norm);
	}
	else if (MethodList.FMH) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing FMH prior gradient");
		*dU += beta * FMH(vec.im_os[0], w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.inffi, w_vec.tr_offsets,
			w_vec.fmh_weights, w_vec.alku_fmh, w_vec.med_no_norm);
	}
	else if (MethodList.WeightedMean) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing weighted mean prior gradient");
		*dU += beta * Weighted_mean(vec.im_os[0], w_vec.Ndx, w_vec.Ndy, w_vec.Ndz, inputScalars, w_vec.weighted_weights,
			w_vec.w_sum, w_vec.mean_type, w_vec.med_no_norm);
	}
	else if (MethodList.TV) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing TV prior gradient");
		status = TVprior(inputScalars, w_vec.data, vec.im_os[0], w_vec, proj, *dU, beta);
	}
	else if (MethodList.hyperbolic) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing hyperbolic prior gradient");
		status = hyperbolic(inputScalars, vec.im_os[0], w_vec, proj, *dU, beta);
	}
	else if (MethodList.AD) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing AD prior gradient");
		if (osa_iter == 0u) {
			*dU += af::constant(0.f, inputScalars.im_dim[0], 1);
		}
		else {
			*dU += AD(vec.im_os[0], inputScalars, w_vec.TimeStepAD, w_vec.KAD, w_vec.NiterAD, w_vec.FluxType,
				w_vec.DiffusionType, w_vec.med_no_norm);
		}
	}
	else if (MethodList.APLS) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing APLS prior gradient");
		status = TVprior(inputScalars, w_vec.data, vec.im_os[0], w_vec, proj, *dU, beta);
	}
	//else if (MethodList.ProxTGV || MethodList.TGV) {
	//	if (inputScalars.verbose >= 3)
	//		mexPrint("Computing TGV prior");
	//	status = TGV(vec.im_os[0], inputScalars, w_vec.data.NiterTGV, w_vec.data.TGVAlpha, w_vec.data.TGVBeta);
	//}
	else if (MethodList.ProxTGV || MethodList.TGV) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing TGV prior");
		status = proxTGV(vec.im_os[0], inputScalars, vec, proj, w_vec, *dU, osa_iter);
	}
	else if (MethodList.ProxTV) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing proximal TV prior");
		status = proxTV(vec.im_os[0], inputScalars, vec, proj, w_vec, *dU, w_vec.betaReg);
	}
	else if (MethodList.NLM) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing NLM prior gradient");
		status = NLM(proj, vec.im_os[0], w_vec, inputScalars, *dU, beta);
	}
	else if (MethodList.RDP) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing RDP prior gradient");
		status = RDP(vec.im_os[0], inputScalars, w_vec.RDP_gamma, proj, *dU, beta, w_vec.RDPref, w_vec.RDPLargeNeighbor, w_vec.RDP_anatomical);
	}
	else if (MethodList.GGMRF) {
		if (inputScalars.verbose >= 3)
			mexPrint("Computing GGMRF prior gradient");
		status = GGMRF(vec.im_os[0], inputScalars, w_vec.GGMRF_p, w_vec.GGMRF_q, w_vec.GGMRF_c, w_vec.GGMRF_pqc, proj, *dU, beta);
	}
	//else if (MethodList.ProxRDP) {
	//	if (inputScalars.verbose >= 3)
	//		mexPrint("Computing proximal RDP prior");
	//	status = proxRDP(vec.im_os[0], inputScalars, vec, proj, *dU, beta, w_vec.RDP_gamma);
	//}
	//else if (MethodList.ProxNLM) {
	//	if (inputScalars.verbose >= 3)
	//		mexPrint("Computing proximal NLM prior");
	//	status = proxNLM(vec.im_os[0], inputScalars, vec, proj, *dU, w_vec, beta);
	//}
	af::deviceGC();
	if (inputScalars.verbose >= 3 && (MethodList.MRP || MethodList.Quad || MethodList.Huber || MethodList.L || MethodList.FMH || MethodList.TV 
		|| MethodList.WeightedMean || MethodList.AD || MethodList.APLS || MethodList.TGV || MethodList.NLM || MethodList.RDP || MethodList.ProxTGV 
		|| MethodList.ProxTV || MethodList.ProxRDP || MethodList.ProxNLM || MethodList.GGMRF))
		mexPrint("Prior computed");
	af::eval(*dU);
	return status;
}